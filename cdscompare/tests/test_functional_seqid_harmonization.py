from pathlib import Path

import pytest

from cdscompare.python_util.locus import (
    STRING_CACHE_DIRECT,
    STRING_CACHE_REVERSE,
    gff_to_cdsInfo,
)
from cdscompare.python_util.seqid_harmonization import (
    attempt_chr_prefix_harmonization,
    check_common_seqids,
)


def write_minimal_gff(path: Path, seqids: list[str]) -> None:
    lines = ["##gff-version 3"]

    for index, seqid in enumerate(seqids, start=1):
        start = index * 1000
        end = start + 299

        gene_id = f"gene_{index}"
        mrna_id = f"mrna_{index}"

        lines.extend(
            [
                f"{seqid}\t.\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}",
                f"{seqid}\t.\tmRNA\t{start}\t{end}\t.\t+\t.\tID={mrna_id};Parent={gene_id}",
                f"{seqid}\t.\tCDS\t{start}\t{end}\t.\t+\t0\tParent={mrna_id}",
            ]
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_and_attempt_harmonization(tmp_path: Path, ref_seqids: list[str], alt_seqids: list[str]):
    ref_gff = tmp_path / "ref.gff3"
    alt_gff = tmp_path / "alt.gff3"

    write_minimal_gff(ref_gff, ref_seqids)
    write_minimal_gff(alt_gff, alt_seqids)

    read_ref = gff_to_cdsInfo(ref_gff)
    read_alt = gff_to_cdsInfo(alt_gff)

    return attempt_chr_prefix_harmonization(read_ref, read_alt)


def seqids_from_read(read: dict[str, list]) -> set[str]:
    seqids = set()

    direct_suffix = "_" + STRING_CACHE_DIRECT
    reverse_suffix = "_" + STRING_CACHE_REVERSE

    for key in read:
        if key.endswith(direct_suffix):
            seqids.add(key[:-len(direct_suffix)])
        elif key.endswith(reverse_suffix):
            seqids.add(key[:-len(reverse_suffix)])
        else:
            raise ValueError(f"Invalid chromosome/strand key in test: {key!r}")

    return seqids


def test_harmonizes_chr_prefix_when_all_seqids_differ(tmp_path, capsys):
    read_ref, read_alt = read_and_attempt_harmonization(
        tmp_path,
        ref_seqids=["chr1", "chr2", "chr3"],
        alt_seqids=["1", "2", "3"],
    )

    assert seqids_from_read(read_ref) == {"1", "2", "3"}
    assert seqids_from_read(read_alt) == {"1", "2", "3"}

    check_common_seqids(read_ref, read_alt)

    captured = capsys.readouterr()
    assert "CDScompare harmonized sequence names" in captured.err
    assert "Shared seqids: 0 before harmonization, 3 after harmonization" in captured.err


def test_harmonizes_chr_prefix_when_seqids_are_mixed(tmp_path, capsys):
    read_ref, read_alt = read_and_attempt_harmonization(
        tmp_path,
        ref_seqids=["chr1", "chr2", "chr3"],
        alt_seqids=["chr1", "2", "chr3"],
    )

    assert seqids_from_read(read_ref) == {"1", "2", "3"}
    assert seqids_from_read(read_alt) == {"1", "2", "3"}

    check_common_seqids(read_ref, read_alt)

    captured = capsys.readouterr()
    assert "CDScompare harmonized sequence names" in captured.err
    assert "Shared seqids: 2 before harmonization, 3 after harmonization" in captured.err


def test_does_not_harmonize_when_it_would_merge_seqids_from_same_annotation(tmp_path, capsys):
    read_ref, read_alt = read_and_attempt_harmonization(
        tmp_path,
        ref_seqids=["chr1", "1", "chr2"],
        alt_seqids=["1", "2", "3"],
    )

    assert seqids_from_read(read_ref) == {"chr1", "1", "chr2"}
    assert seqids_from_read(read_alt) == {"1", "2", "3"}

    check_common_seqids(read_ref, read_alt)

    captured = capsys.readouterr()
    assert captured.err == ""


def test_does_not_warn_when_seqids_already_match(tmp_path, capsys):
    read_ref, read_alt = read_and_attempt_harmonization(
        tmp_path,
        ref_seqids=["chr1", "chr2"],
        alt_seqids=["chr1", "chr2"],
    )

    assert seqids_from_read(read_ref) == {"chr1", "chr2"}
    assert seqids_from_read(read_alt) == {"chr1", "chr2"}

    check_common_seqids(read_ref, read_alt)

    captured = capsys.readouterr()
    assert captured.err == ""


def test_raises_explicit_error_when_no_common_seqid_after_harmonization_attempt(tmp_path, capsys):
    read_ref, read_alt = read_and_attempt_harmonization(
        tmp_path,
        ref_seqids=["chr1", "chr2"],
        alt_seqids=["scaffoldA", "scaffoldB"],
    )

    with pytest.raises(SystemExit, match="No common seqids were found between annotations"):
        check_common_seqids(read_ref, read_alt)

    captured = capsys.readouterr()
    assert captured.err == ""
