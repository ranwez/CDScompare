# -*- coding: utf-8 -*-

from collections import defaultdict
import sys

from cdscompare.python_util.locus import STRING_CACHE_DIRECT, STRING_CACHE_REVERSE


def _split_chr_strand_key(key: str) -> tuple[str, str]:
    direct_suffix = "_" + STRING_CACHE_DIRECT
    reverse_suffix = "_" + STRING_CACHE_REVERSE

    if key.endswith(direct_suffix):
        return key[:-len(direct_suffix)], STRING_CACHE_DIRECT

    if key.endswith(reverse_suffix):
        return key[:-len(reverse_suffix)], STRING_CACHE_REVERSE

    raise ValueError(f"Invalid chromosome/strand key: {key!r}")


def _harmonize_seqid(seqid: str) -> str:
    if seqid[:3].lower() == "chr":
        return seqid[3:]

    return seqid


def attempt_chr_prefix_harmonization(
    read_ref: dict[str, list],
    read_alt: dict[str, list],
) -> tuple[dict[str, list], dict[str, list]]:
    ref_seqids = {_split_chr_strand_key(key)[0] for key in read_ref}
    alt_seqids = {_split_chr_strand_key(key)[0] for key in read_alt}

    common_before = len(ref_seqids & alt_seqids)

    ref_seqid_map = {seqid: _harmonize_seqid(seqid) for seqid in ref_seqids}
    alt_seqid_map = {seqid: _harmonize_seqid(seqid) for seqid in alt_seqids}

    # If several seqids from the same annotation would become identical,
    # harmonization is ambiguous. In that case, keep original names.
    if len(set(ref_seqid_map.values())) != len(ref_seqids) or len(set(alt_seqid_map.values())) != len(alt_seqids):
        return read_ref, read_alt

    candidate_ref_seqids = set(ref_seqid_map.values())
    candidate_alt_seqids = set(alt_seqid_map.values())

    common_after = len(candidate_ref_seqids & candidate_alt_seqids)

    if common_after <= common_before:
        return read_ref, read_alt

    harmonized_ref = defaultdict(list)
    harmonized_alt = defaultdict(list)

    for key, loci in read_ref.items():
        seqid, direction = _split_chr_strand_key(key)
        harmonized_ref[f"{ref_seqid_map[seqid]}_{direction}"] = loci

    for key, loci in read_alt.items():
        seqid, direction = _split_chr_strand_key(key)
        harmonized_alt[f"{alt_seqid_map[seqid]}_{direction}"] = loci

    print(
        "[WARNING] CDScompare harmonized sequence names by removing 'chr' prefixes. "
        f"Shared seqids: {common_before} before harmonization, "
        f"{common_after} after harmonization.",
        file=sys.stderr,
    )

    return harmonized_ref, harmonized_alt



def check_common_seqids(
    read_ref: dict[str, list],
    read_alt: dict[str, list],
) -> None:
    ref_seqids = {_split_chr_strand_key(key)[0] for key in read_ref}
    alt_seqids = {_split_chr_strand_key(key)[0] for key in read_alt}

    common = ref_seqids & alt_seqids

    if common:
        return

    raise SystemExit(
        "\nNo common seqids were found between annotations. "
        "Please check that both GFF files use compatible seqids in column 1. "
        f"\nReference seqid examples: {sorted(ref_seqids)[:5]}; "
        f"\nAlternative seqid examples: {sorted(alt_seqids)[:5]}"
    )
