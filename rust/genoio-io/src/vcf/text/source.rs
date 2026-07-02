//! Reader and index setup for the text VCF backend.
//!
//! This module is the imperative shell for text VCF reads: it opens files,
//! configures optional BGZF threading, reads the sample header line, and maps
//! concrete region predicates to tabix/CSI chunks.

// pattern: Imperative Shell

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::NonZero;
use std::path::Path;

use flate2::read::MultiGzDecoder;
use genoio_core::{
    select_samples_source_order, DenseSampleSelection, GenoioError, RegionPredicate,
};
use noodles_bgzf as bgzf;
use noodles_core::{Position, Region};
use noodles_csi::{self as csi, BinningIndex};
use noodles_tabix as tabix;
use noodles_vcf as noodles;

use crate::error::Result;

use super::header::read_sample_records_from_header;
use super::VCF_TEXT_BUFFER_SIZE;
use crate::vcf::{is_bcf_path, is_compressed_vcf, policy::companion_index_path};

pub(super) type CompressedVcfReader = noodles::io::Reader<BufReader<MultiGzDecoder<File>>>;
pub(super) type ThreadedCompressedVcfReader =
    noodles::io::Reader<bgzf::io::MultithreadedReader<File>>;
type IndexedCompressedVcfReader<'a, R> = noodles::io::Reader<csi::io::Query<'a, R>>;
pub(super) type IndexedBgzfReader<'a> = IndexedCompressedVcfReader<'a, bgzf::io::Reader<File>>;
pub(super) type ThreadedIndexedBgzfReader<'a> =
    IndexedCompressedVcfReader<'a, bgzf::io::MultithreadedReader<File>>;
type IndexChunk = csi::binning_index::index::reference_sequence::bin::Chunk;
pub(super) type PlainVcfReader = noodles::io::Reader<BufReader<File>>;

pub(super) struct TextVcfInput<R> {
    pub(super) reader: noodles::io::Reader<R>,
    pub(super) selection: DenseSampleSelection,
}

// Dispatch once at setup so hot record loops stay monomorphized for the
// underlying reader type instead of paying through `dyn BufRead`.
pub(super) enum TextVcfSource {
    Compressed(TextVcfInput<BufReader<MultiGzDecoder<File>>>),
    ThreadedCompressed(TextVcfInput<bgzf::io::MultithreadedReader<File>>),
    Plain(TextVcfInput<BufReader<File>>),
}

#[derive(Clone, Copy)]
pub(super) struct DenseReadSource<'a> {
    pub(super) region: Option<&'a RegionPredicate>,
}

impl<'a> DenseReadSource<'a> {
    pub(super) const fn full_scan() -> Self {
        Self { region: None }
    }
}

pub(super) struct IndexedTextVcfInput<'a> {
    pub(super) selection: DenseSampleSelection,
    pub(super) region: &'a RegionPredicate,
}

impl<'a> IndexedTextVcfInput<'a> {
    pub(super) const fn dense_source(&self) -> DenseReadSource<'a> {
        DenseReadSource {
            region: Some(self.region),
        }
    }
}

pub(super) fn ensure_text_vcf_supported(path: &Path, threads: Option<usize>) -> Result<()> {
    if !is_bcf_path(path) && text_backend_supports_threads(path, threads) {
        return Ok(());
    }
    Err(GenoioError::internal_contract(format!(
        "text VCF backend called for unsupported source {}",
        path.display()
    )))
}

pub(super) fn ensure_text_indexed_vcf_supported(path: &Path) -> Result<()> {
    if is_compressed_vcf(path) {
        return Ok(());
    }
    Err(GenoioError::internal_contract(format!(
        "indexed text VCF backend called for uncompressed source {}",
        path.display()
    )))
}

fn text_backend_supports_threads(path: &Path, threads: Option<usize>) -> bool {
    // Noodles only gives us threaded BGZF decompression. Plain text with an
    // explicit thread count is rejected at the public VCF boundary.
    threads.is_none() || is_compressed_vcf(path)
}

pub(super) fn open_text_vcf_input(
    path: &Path,
    requested_samples: Option<&[String]>,
    threads: Option<usize>,
) -> Result<TextVcfSource> {
    if is_compressed_vcf(path) {
        match threads {
            Some(threads) => open_text_vcf_input_from_reader(
                path,
                requested_samples,
                open_threaded_compressed_reader(path, threads)?,
            )
            .map(TextVcfSource::ThreadedCompressed),
            None => open_text_vcf_input_from_reader(
                path,
                requested_samples,
                open_compressed_reader(path)?,
            )
            .map(TextVcfSource::Compressed),
        }
    } else {
        open_text_vcf_input_from_reader(path, requested_samples, open_plain_reader(path)?)
            .map(TextVcfSource::Plain)
    }
}

fn open_text_vcf_input_from_reader<R: BufRead>(
    path: &Path,
    requested_samples: Option<&[String]>,
    mut reader: noodles::io::Reader<R>,
) -> Result<TextVcfInput<R>> {
    let all_samples = read_sample_records_from_header(path, reader.get_mut())?;
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    if requested_samples.is_some() {
        for (sample, source_index) in selection.samples.iter_mut().zip(&selection.source_indices) {
            sample.source_sample_index = Some(*source_index);
        }
    }
    Ok(TextVcfInput { reader, selection })
}

pub(super) fn open_text_sample_selection(
    path: &Path,
    requested_samples: Option<&[String]>,
    threads: Option<usize>,
) -> Result<DenseSampleSelection> {
    Ok(open_text_vcf_input(path, requested_samples, threads)?.into_selection())
}

pub(super) fn with_indexed_text_vcf_input<'region, T, ReadRecords, EmptyResult>(
    path: &Path,
    requested_samples: Option<&[String]>,
    region: &'region RegionPredicate,
    read_records: ReadRecords,
    empty_result: EmptyResult,
) -> Result<T>
where
    ReadRecords: for<'reader> FnOnce(
        IndexedTextVcfInput<'region>,
        &mut IndexedBgzfReader<'reader>,
    ) -> Result<T>,
    EmptyResult: FnOnce(DenseSampleSelection) -> Result<T>,
{
    run_indexed_text_vcf_input(
        path,
        requested_samples,
        region,
        open_bgzf_reader(path)?,
        read_records,
        empty_result,
    )
}

pub(super) fn with_threaded_indexed_text_vcf_input<'region, T, ReadRecords, EmptyResult>(
    path: &Path,
    requested_samples: Option<&[String]>,
    region: &'region RegionPredicate,
    threads: usize,
    read_records: ReadRecords,
    empty_result: EmptyResult,
) -> Result<T>
where
    ReadRecords: for<'reader> FnOnce(
        IndexedTextVcfInput<'region>,
        &mut ThreadedIndexedBgzfReader<'reader>,
    ) -> Result<T>,
    EmptyResult: FnOnce(DenseSampleSelection) -> Result<T>,
{
    run_indexed_text_vcf_input(
        path,
        requested_samples,
        region,
        open_threaded_bgzf_reader(path, threads)?,
        read_records,
        empty_result,
    )
}

fn run_indexed_text_vcf_input<'region, T, R, ReadRecords, EmptyResult>(
    path: &Path,
    requested_samples: Option<&[String]>,
    region: &'region RegionPredicate,
    mut bgzf_reader: R,
    read_records: ReadRecords,
    empty_result: EmptyResult,
) -> Result<T>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
    ReadRecords: for<'reader> FnOnce(
        IndexedTextVcfInput<'region>,
        &mut IndexedCompressedVcfReader<'reader, R>,
    ) -> Result<T>,
    EmptyResult: FnOnce(DenseSampleSelection) -> Result<T>,
{
    let chunks = index_chunks_for_region(path, region)?;
    let all_samples = read_sample_records_from_header(path, &mut bgzf_reader)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;

    let Some(chunks) = chunks else {
        return empty_result(selection);
    };

    // Query chunks directly rather than using noodles' strict IndexedReader.
    // This path only needs the #CHROM line, and real VCFs can have header
    // records that permissive parsers accept but noodles' full header parser
    // rejects.
    let query = csi::io::Query::new(&mut bgzf_reader, chunks);
    let mut reader = noodles::io::Reader::new(query);
    read_records(IndexedTextVcfInput { selection, region }, &mut reader)
}

fn open_bgzf_reader(path: &Path) -> Result<bgzf::io::Reader<File>> {
    File::open(path)
        .map(bgzf::io::Reader::new)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))
}

fn open_threaded_bgzf_reader(
    path: &Path,
    threads: usize,
) -> Result<bgzf::io::MultithreadedReader<File>> {
    let worker_count = NonZero::new(threads).ok_or_else(|| {
        GenoioError::invalid_source(path, "vcf thread count must be greater than zero")
    })?;
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))?;
    // Thread only BGZF inflation. Record parsing stays ordered and uses the
    // same decode code as the unthreaded text backend.
    Ok(bgzf::io::MultithreadedReader::with_worker_count(
        worker_count,
        file,
    ))
}

fn index_chunks_for_region(
    path: &Path,
    region: &RegionPredicate,
) -> Result<Option<Vec<IndexChunk>>> {
    let index = read_associated_index(path).map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf index read error: {error}"))
    })?;
    let region = noodles_region_from_predicate(path, region)?;
    let Some(header) = index.header() else {
        return Ok(None);
    };
    let Some(reference_sequence_id) = header
        .reference_sequence_names()
        .get_index_of(region.name())
    else {
        return Ok(None);
    };
    index
        .query(reference_sequence_id, region.interval())
        .map(Some)
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf index query error: {error}"))
        })
}

fn read_associated_index(path: &Path) -> std::io::Result<Box<dyn BinningIndex>> {
    match tabix::fs::read(companion_index_path(path, "tbi")) {
        Ok(index) => Ok(Box::new(index)),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            csi::fs::read(companion_index_path(path, "csi"))
                .map(|index| Box::new(index) as Box<dyn BinningIndex>)
        }
        Err(error) => Err(error),
    }
}

fn noodles_region_from_predicate(path: &Path, region: &RegionPredicate) -> Result<Region> {
    let start = position_from_u32(path, region.start, "start")?;
    let end = position_from_u32(path, region.end, "end")?;
    Ok(Region::new(region.chrom.as_str(), start..=end))
}

fn position_from_u32(path: &Path, value: u32, label: &str) -> Result<Position> {
    let value = usize::try_from(value).map_err(|_| {
        GenoioError::invalid_source(
            path,
            format!("vcf region {label} coordinate is out of range"),
        )
    })?;
    Position::try_from(value).map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("vcf region {label} coordinate is invalid: {error}"),
        )
    })
}

pub(super) fn open_compressed_reader(path: &Path) -> Result<CompressedVcfReader> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))?;
    let reader = BufReader::with_capacity(VCF_TEXT_BUFFER_SIZE, MultiGzDecoder::new(file));
    Ok(noodles::io::Reader::new(reader))
}

fn open_threaded_compressed_reader(
    path: &Path,
    threads: usize,
) -> Result<ThreadedCompressedVcfReader> {
    open_threaded_bgzf_reader(path, threads).map(noodles::io::Reader::new)
}

pub(super) fn open_plain_reader(path: &Path) -> Result<PlainVcfReader> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))?;
    let reader = BufReader::with_capacity(VCF_TEXT_BUFFER_SIZE, file);
    Ok(noodles::io::Reader::new(reader))
}
