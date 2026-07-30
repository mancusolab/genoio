# pattern: Imperative Shell

from __future__ import annotations

import tempfile
from importlib.resources import files
from pathlib import Path

import numpy as np

import genoio

EXPECTED_MATRIX = np.array(
    [
        [0.0, np.nan, 2.0, 0.0],
        [np.nan, 0.0, 1.0, 1.0],
        [2.0, 1.0, 0.0, 2.0],
        [1.0, 2.0, np.nan, 0.0],
    ],
    dtype=np.float32,
)


def _assert_shape(name: str, actual: tuple[int, ...], expected: tuple[int, ...]) -> None:
    if actual != expected:
        raise AssertionError(f"{name} shape mismatch: expected {expected}, got {actual}")


def _assert_close(name: str, actual: np.ndarray, expected: np.ndarray) -> None:
    if actual.shape != expected.shape:
        raise AssertionError(f"{name} shape mismatch: expected {expected.shape}, got {actual.shape}")
    np.testing.assert_allclose(actual, expected, equal_nan=True)


def _smoke_vcf(tmp_path: Path) -> None:
    matrix, variants = genoio.vcf(_write_vcf(tmp_path)).read(missing="nan", return_variants=True)
    _assert_close("vcf", matrix, EXPECTED_MATRIX)
    if variants["id"].to_list() != ["rs1", "rs2", "indel1", "rs4"]:
        raise AssertionError("vcf variant ids do not match expected tiny fixture")


def _smoke_plink1(tmp_path: Path) -> None:
    matrix, samples = genoio.bfile(_write_plink1(tmp_path)).read(missing="nan", return_samples=True)
    _assert_close("plink1", matrix, EXPECTED_MATRIX)
    if samples["iid"].to_list() != ["S1", "S2", "S3", "S4"]:
        raise AssertionError("plink1 sample ids do not match expected tiny fixture")


def _smoke_plink2(tmp_path: Path) -> None:
    matrix = genoio.pfile(_write_plink2(tmp_path)).read(missing="nan")
    _assert_close("plink2", matrix, EXPECTED_MATRIX)


def _smoke_bgen(tmp_path: Path) -> None:
    matrix, variants = genoio.bgen(_write_bgen(tmp_path)).read(
        dosage="dosage",
        missing="nan",
        return_variants=True,
    )
    _assert_shape("bgen", matrix.shape, (2, 2))
    if variants["id"].to_list() != ["rs1", "rs2"]:
        raise AssertionError("bgen variant ids do not match expected tiny fixture")


def main() -> None:
    print(f"imported genoio from {genoio.__file__}")
    if not files("genoio").joinpath("py.typed").is_file():
        raise AssertionError("installed genoio wheel is missing py.typed")
    with tempfile.TemporaryDirectory(prefix="genoio-wheel-smoke-") as tmpdir:
        tmp_path = Path(tmpdir)
        _smoke_vcf(tmp_path)
        _smoke_plink1(tmp_path)
        _smoke_plink2(tmp_path)
        _smoke_bgen(tmp_path)
    print("installed wheel IO smoke passed")


def _write_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "canonical.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
1\t10\trs1\tG\tA\t.\tPASS\t.\tGT\t0/0\t./.\t1/1\t0/1
1\t20\trs2\tT\tC\t.\tPASS\t.\tGT\t./.\t0/0\t0/1\t1/1
2\t30\tindel1\tA\tAT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\t./.
2\t40\trs4\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
"""
    )
    return path


def _write_plink1(tmp_path: Path) -> Path:
    prefix = tmp_path / "canonical_bed"
    prefix.with_suffix(".bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x87, 0x2D, 0x78, 0xCB]))
    prefix.with_suffix(".bim").write_text(
        """\
1 rs1 0 10 A G
1 rs2 0 20 C T
2 indel1 0 30 AT A
2 rs4 0 40 T C
"""
    )
    prefix.with_suffix(".fam").write_text(
        """\
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
F2 S4 0 0 2 -9
"""
    )
    return prefix


def _write_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "canonical_pgen"
    prefix.with_suffix(".pgen").write_bytes(
        bytes([0x6C, 0x1B, 0x02])
        + (4).to_bytes(4, "little")
        + (4).to_bytes(4, "little")
        + bytes([0x00, 0x6C, 0x93, 0xC6, 0x24])
    )
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT
1 10 rs1 G A
1 20 rs2 T C
2 30 indel1 A AT
2 40 rs4 C T
"""
    )
    prefix.with_suffix(".psam").write_text(
        """\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
F2 S4 0 0 2 -9
"""
    )
    return prefix


def _write_bgen(tmp_path: Path) -> Path:
    path = tmp_path / "dosage.bgen"
    contents = bytearray()
    variants = [
        ("var1", "rs1", "1", 10, ["A", "G"]),
        ("var2", "rs2", "2", 20, ["C", "T"]),
    ]
    variant_calls = [
        [(204, 26), (51, 128)],
        [(0, 255), (102, 102)],
    ]
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((20).to_bytes(4, "little"))
    contents.extend(len(variants).to_bytes(4, "little"))
    contents.extend((2).to_bytes(4, "little"))
    contents.extend(b"bgen")
    contents.extend(((2 << 2) | (1 << 31)).to_bytes(4, "little"))
    contents.extend(_bgen_sample_identifier_block(["sample_1", "sample_2"]))
    contents[0:4] = (len(contents) - 4).to_bytes(4, "little")

    for variant, calls in zip(variants, variant_calls, strict=True):
        contents.extend(_bgen_variant_identifying_data(*variant))
        contents.extend(_bgen_dosage_probability_block(calls))

    path.write_bytes(contents)
    return path


def _bgen_sample_identifier_block(sample_ids: list[str]) -> bytes:
    contents = bytearray()
    block_len = 8 + sum(2 + len(sample_id.encode()) for sample_id in sample_ids)
    contents.extend(block_len.to_bytes(4, "little"))
    contents.extend(len(sample_ids).to_bytes(4, "little"))
    for sample_id in sample_ids:
        encoded = sample_id.encode()
        contents.extend(len(encoded).to_bytes(2, "little"))
        contents.extend(encoded)
    return bytes(contents)


def _bgen_variant_identifying_data(
    variant_id: str,
    rsid: str,
    chrom: str,
    pos: int,
    alleles: list[str],
) -> bytes:
    contents = bytearray()
    for value in (variant_id, rsid, chrom):
        encoded = value.encode()
        contents.extend(len(encoded).to_bytes(2, "little"))
        contents.extend(encoded)
    contents.extend(pos.to_bytes(4, "little"))
    contents.extend(len(alleles).to_bytes(2, "little"))
    for allele in alleles:
        encoded = allele.encode()
        contents.extend(len(encoded).to_bytes(4, "little"))
        contents.extend(encoded)
    return bytes(contents)


def _bgen_dosage_probability_block(calls: list[tuple[int, int]]) -> bytes:
    payload = bytearray()
    payload.extend(len(calls).to_bytes(4, "little"))
    payload.extend((2).to_bytes(2, "little"))
    payload.extend((2).to_bytes(1, "little"))
    payload.extend((2).to_bytes(1, "little"))
    payload.extend(2 for _ in calls)
    payload.extend((0).to_bytes(1, "little"))
    payload.extend((8).to_bytes(1, "little"))
    _append_bgen_packed_probabilities(payload, calls)

    contents = bytearray()
    contents.extend(len(payload).to_bytes(4, "little"))
    contents.extend(payload)
    return bytes(contents)


def _append_bgen_packed_probabilities(output: bytearray, calls: list[tuple[int, int]]) -> None:
    current_byte = 0
    bits_in_current_byte = 0
    for call in calls:
        for value in call:
            for bit_index in range(8):
                current_byte |= ((value >> bit_index) & 1) << bits_in_current_byte
                bits_in_current_byte += 1
                if bits_in_current_byte == 8:
                    output.append(current_byte)
                    current_byte = 0
                    bits_in_current_byte = 0
    if bits_in_current_byte:
        output.append(current_byte)


if __name__ == "__main__":
    main()
