# pattern: Imperative Shell

from pathlib import Path

import numpy as np

EXPECTED_MATRIX = np.array(
    [
        [0.0, np.nan, 2.0, 0.0],
        [np.nan, 0.0, 1.0, 1.0],
        [2.0, 1.0, 0.0, 2.0],
        [1.0, 2.0, np.nan, 0.0],
    ],
    dtype=np.float32,
)
EXPECTED_SAMPLES = ["S1", "S2", "S3", "S4"]
EXPECTED_VARIANTS = ["rs1", "rs2", "indel1", "rs4"]
EXPECTED_VARIANT_ROWS = [
    ("1", 10, "rs1", "G", "A"),
    ("1", 20, "rs2", "T", "C"),
    ("2", 30, "indel1", "A", "AT"),
    ("2", 40, "rs4", "C", "T"),
]


def write_bgen_dosage(
    tmp_path: Path,
    *,
    missing: bool = False,
    phased: bool = False,
    pack_missing_probabilities: bool = False,
    sample_ids: list[str] | None = None,
    variant_calls: list[list[tuple[int, int] | None]] | None = None,
    variants: list[tuple[str, str, str, int, list[str]]] | None = None,
) -> Path:
    """Write a tiny Layout 2 BGEN dosage fixture.

    `pack_missing_probabilities` mirrors writers that emit zero-valued packed
    probabilities for samples already marked missing in the ploidy bytes.
    """
    path = tmp_path / "dosage.bgen"
    contents = bytearray()
    flags = (2 << 2) | (1 << 31)
    sample_ids = ["sample_1", "sample_2"] if sample_ids is None else sample_ids
    variants = (
        [
            ("var1", "rs1", "1", 10, ["A", "G"]),
            ("var2", "rs2", "2", 20, ["C", "T"]),
        ]
        if variants is None
        else variants
    )
    if variant_calls is None:
        variant_calls = [
            [(204, 26), None if missing else (51, 128)],
            [(0, 255), (102, 102)],
        ]
    if len(variant_calls) != len(variants):
        raise ValueError("variant_calls and variants must have matching lengths")

    contents.extend((20).to_bytes(4, "little"))
    contents.extend((20).to_bytes(4, "little"))
    contents.extend(len(variants).to_bytes(4, "little"))
    contents.extend(len(sample_ids).to_bytes(4, "little"))
    contents.extend(b"bgen")
    contents.extend(flags.to_bytes(4, "little"))
    contents.extend(_bgen_sample_identifier_block(sample_ids))
    variant_offset = len(contents) - 4
    contents[0:4] = variant_offset.to_bytes(4, "little")

    for variant, calls in zip(variants, variant_calls, strict=True):
        contents.extend(_bgen_variant_identifying_data(*variant))
        contents.extend(
            _bgen_dosage_probability_block(
                8,
                calls,
                phased=phased,
                pack_missing_probabilities=pack_missing_probabilities,
            )
        )

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


def _bgen_dosage_probability_block(
    bit_depth: int,
    calls: list[tuple[int, int] | None],
    *,
    phased: bool = False,
    pack_missing_probabilities: bool = False,
) -> bytes:
    payload = bytearray()
    payload.extend(len(calls).to_bytes(4, "little"))
    payload.extend((2).to_bytes(2, "little"))
    payload.extend((2).to_bytes(1, "little"))
    payload.extend((2).to_bytes(1, "little"))
    payload.extend((2 if call is not None else 0b1000_0010) for call in calls)
    payload.extend((1 if phased else 0).to_bytes(1, "little"))
    payload.extend(bit_depth.to_bytes(1, "little"))
    packed_calls = _packed_bgen_calls(calls, pack_missing_probabilities=pack_missing_probabilities)
    _append_bgen_packed_probabilities(payload, bit_depth, packed_calls)

    contents = bytearray()
    contents.extend(len(payload).to_bytes(4, "little"))
    contents.extend(payload)
    return bytes(contents)


def _packed_bgen_calls(
    calls: list[tuple[int, int] | None],
    *,
    pack_missing_probabilities: bool,
) -> list[tuple[int, int]]:
    if pack_missing_probabilities:
        return [(0, 0) if call is None else call for call in calls]
    return [call for call in calls if call is not None]


def _append_bgen_packed_probabilities(
    output: bytearray,
    bit_depth: int,
    calls: list[tuple[int, int]],
) -> None:
    current_byte = 0
    bits_in_current_byte = 0
    for call in calls:
        for value in call:
            for bit_index in range(bit_depth):
                current_byte |= ((value >> bit_index) & 1) << bits_in_current_byte
                bits_in_current_byte += 1
                if bits_in_current_byte == 8:
                    output.append(current_byte)
                    current_byte = 0
                    bits_in_current_byte = 0
    if bits_in_current_byte:
        output.append(current_byte)


def write_fixed_width_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "tiny"
    prefix.with_suffix(".pgen").write_bytes(
        bytes(
            [
                0x6C,
                0x1B,
                0x02,
                0x03,
                0x00,
                0x00,
                0x00,
                0x03,
                0x00,
                0x00,
                0x00,
                0x00,
                0x2C,
                0x11,
                0x06,
            ]
        )
    )
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
2 30 rs3 G A 50
"""
    )
    prefix.with_suffix(".psam").write_text(
        """\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
"""
    )
    return prefix


def write_fixed_width_plink2_dosage(tmp_path: Path) -> Path:
    prefix = tmp_path / "dosage"

    def scaled(value: float) -> bytes:
        return round(value / 2.0 * 32768.0).to_bytes(2, "little")

    variant_1_hardcalls = bytes([0x24])
    variant_2_hardcalls = bytes([0x0C])
    prefix.with_suffix(".pgen").write_bytes(
        bytes([0x6C, 0x1B, 0x03])
        + (2).to_bytes(4, "little")
        + (3).to_bytes(4, "little")
        + bytes([0])
        + variant_1_hardcalls
        + scaled(0.2)
        + scaled(1.4)
        + scaled(1.8)
        + variant_2_hardcalls
        + scaled(0.0)
        + (65535).to_bytes(2, "little")
        + scaled(0.7)
    )
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
"""
    )
    prefix.with_suffix(".psam").write_text(
        """\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
"""
    )
    return prefix


def write_fixed_width_phased_dosage_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "fixed_width_phased_dosage"
    record_1 = (
        bytes([0x25])
        + _plink2_scaled_dosage(1.0)
        + _plink2_scaled_dosage(0.5)
        + _plink2_scaled_dosage(2.0)
        + _plink2_scaled_phase_delta(0.25, 0.75)
        + _plink2_scaled_phase_delta(0.0, 0.5)
        + _plink2_scaled_phase_delta(1.0, 1.0)
    )
    record_2 = (
        bytes([0x00])
        + _plink2_scaled_dosage(0.0)
        + _plink2_scaled_dosage(0.2)
        + _plink2_scaled_dosage(0.4)
        + _plink2_scaled_phase_delta(0.0, 0.0)
        + _plink2_scaled_phase_delta(0.1, 0.1)
        + _plink2_scaled_phase_delta(0.2, 0.2)
    )
    prefix.with_suffix(".pgen").write_bytes(
        bytes([0x6C, 0x1B, 0x04])
        + (2).to_bytes(4, "little")
        + (3).to_bytes(4, "little")
        + bytes([0])
        + record_1
        + record_2
    )
    _write_plink2_pvar(prefix)
    _write_plink2_psam(prefix)
    return prefix


def write_phased_hardcall_plink2(tmp_path: Path, *, unphased_second_variant: bool = False) -> Path:
    prefix = tmp_path / "phased_hardcall"
    record_1 = bytes([0x21, 0x00])
    record_2 = bytes([0x35, 0x02])
    record_types = [0x10, 0x00 if unphased_second_variant else 0x10]
    records = [record_1, record_2[:-1] if unphased_second_variant else record_2]
    _write_variable_width_plink2(prefix, record_types, records, n_samples=3)
    _write_plink2_pvar(prefix)
    _write_plink2_psam(prefix)
    return prefix


def write_phased_dosage_plink2(tmp_path: Path, *, unphased_second_variant: bool = False) -> Path:
    prefix = tmp_path / "phased_dosage"
    record_1 = (
        bytes([0x25])
        + _plink2_scaled_dosage(1.0)
        + _plink2_scaled_dosage(0.5)
        + _plink2_scaled_dosage(2.0)
        + _plink2_scaled_phase_delta(0.25, 0.75)
        + _plink2_scaled_phase_delta(0.0, 0.5)
        + _plink2_scaled_phase_delta(1.0, 1.0)
    )
    record_2 = (
        bytes([0x00])
        + _plink2_scaled_dosage(0.0)
        + _plink2_scaled_dosage(0.2)
        + _plink2_scaled_dosage(0.4)
        + _plink2_scaled_phase_delta(0.0, 0.0)
        + _plink2_scaled_phase_delta(0.1, 0.1)
        + _plink2_scaled_phase_delta(0.2, 0.2)
    )
    if unphased_second_variant:
        record_2 = bytes([0x00]) + _plink2_scaled_dosage(0.0) + _plink2_scaled_dosage(0.2) + _plink2_scaled_dosage(0.4)
    _write_variable_width_plink2(
        prefix,
        [0xC0, 0x40 if unphased_second_variant else 0xC0],
        [record_1, record_2],
        n_samples=3,
    )
    _write_plink2_pvar(prefix)
    _write_plink2_psam(prefix)
    return prefix


def write_ld_phased_hardcall_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "ld_phased_hardcall"
    record_1 = bytes([0x21, 0x00])
    record_2 = bytes([0x02, 0x01, 0x0D, 0x01, 0x02])
    _write_variable_width_plink2(prefix, [0x10, 0x12], [record_1, record_2], n_samples=3)
    _write_plink2_pvar(prefix)
    _write_plink2_psam(prefix)
    return prefix


def write_ld_phased_dosage_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "ld_phased_dosage"
    record_1 = (
        bytes([0x25])
        + _plink2_scaled_dosage(1.0)
        + _plink2_scaled_dosage(0.5)
        + _plink2_scaled_dosage(2.0)
        + _plink2_scaled_phase_delta(0.25, 0.75)
        + _plink2_scaled_phase_delta(0.0, 0.5)
        + _plink2_scaled_phase_delta(1.0, 1.0)
    )
    record_2 = (
        bytes([0x03, 0x00, 0x00, 0x01, 0x01])
        + _plink2_scaled_dosage(0.0)
        + _plink2_scaled_dosage(0.2)
        + _plink2_scaled_dosage(0.4)
        + _plink2_scaled_phase_delta(0.0, 0.0)
        + _plink2_scaled_phase_delta(0.1, 0.1)
        + _plink2_scaled_phase_delta(0.2, 0.2)
    )
    _write_variable_width_plink2(prefix, [0xC0, 0xC2], [record_1, record_2], n_samples=3)
    _write_plink2_pvar(prefix)
    _write_plink2_psam(prefix)
    return prefix


def write_sample_filtered_unphased_hardcall_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "sample_filtered_unphased_hardcall"
    _write_variable_width_plink2(prefix, [0x10], [bytes([0x15, 0x0D, 0x02])], n_samples=3)
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
"""
    )
    _write_plink2_psam(prefix)
    return prefix


def _write_variable_width_plink2(
    prefix: Path,
    record_types: list[int],
    records: list[bytes],
    *,
    n_samples: int,
) -> None:
    header_len = 12 + 8 + len(record_types) + len(records)
    payload = bytearray([0x6C, 0x1B, 0x10])
    payload.extend(len(records).to_bytes(4, "little"))
    payload.extend(n_samples.to_bytes(4, "little"))
    payload.append(0x04)
    payload.extend(header_len.to_bytes(8, "little"))
    payload.extend(record_types)
    payload.extend(len(record) for record in records)
    for record in records:
        payload.extend(record)
    prefix.with_suffix(".pgen").write_bytes(bytes(payload))


def _write_plink2_pvar(prefix: Path) -> None:
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
"""
    )


def _write_plink2_psam(prefix: Path) -> None:
    prefix.with_suffix(".psam").write_text(
        """\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
"""
    )


def _plink2_scaled_dosage(value: float) -> bytes:
    return round(value / 2.0 * 32768.0).to_bytes(2, "little")


def _plink2_scaled_phase_delta(left: float, right: float) -> bytes:
    raw = round((left - right) * 16384.0)
    return raw.to_bytes(2, "little", signed=True)


def write_canonical_vcf(tmp_path: Path) -> Path:
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


def write_canonical_plink1(tmp_path: Path) -> Path:
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


def write_canonical_plink2(tmp_path: Path) -> Path:
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
