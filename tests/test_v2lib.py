"""Tests for the results/v2 summary logic.

Hand-built inputs with known answers. The read classifier carries the most
weight: it is what the aligner-disagreement explanation rests on, so each
category gets its own case, including the ones that are easy to confuse.
"""
import json

import pytest

import v2lib

FLAGSTAT = """\
120 + 0 in total (QC-passed reads + QC-failed reads)
118 + 0 primary
2 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
120 + 0 mapped (100.00% : N/A)
118 + 0 paired in sequencing
59 + 0 read1
59 + 0 read2
117 + 0 properly paired (99.15% : N/A)
"""


def sam(qname="r1", flag=0x1 | 0x2 | 0x40, rname="NC_1", pos=100, mapq=42,
        cigar="150M", nm=0):
    return v2lib.parse_sam_line(
        f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t=\t300\t350\tACGT\tIIII\tNM:i:{nm}"
    )


# ---------------------------------------------------------------- parsers

def test_flagstat_pairs_halves_and_floors():
    # 117 properly paired reads -> 58 pairs, same floor rule as v1
    assert v2lib.flagstat_proper_pairs(FLAGSTAT) == 58


def test_flagstat_empty_bam_is_zero():
    assert v2lib.flagstat_proper_pairs("0 + 0 in total (QC-passed)\n0 + 0 properly paired (N/A : N/A)\n") == 0


def test_fastp_pairs():
    d = {"summary": {"after_filtering": {"total_reads": 1911862}}}
    assert v2lib.fastp_post_qc_pairs(json.dumps(d)) == 955931


def test_hostile_stats_list_form():
    d = [{"reads_in": 200, "reads_out": 150, "reads_removed": 50,
          "reads_removed_proportion": 0.25}]
    assert v2lib.hostile_stats(json.dumps(d)) == {
        "reads_in": 200, "reads_out": 150, "reads_removed": 50}


def test_hostile_stats_rejects_multiple_records():
    d = [{"reads_in": 1, "reads_out": 1, "reads_removed": 0}] * 2
    with pytest.raises(ValueError):
        v2lib.hostile_stats(json.dumps(d))


def test_parse_sam_soft_clip_and_mate():
    r = sam(flag=0x1 | 0x2 | 0x80, cigar="10S130M10S", nm=3)
    assert r["soft_clipped"] == 20
    assert r["mate"] == 2
    assert r["nm"] == 3


# ---------------------------------------------------------------- host depletion

def _long(rows):
    return {(s, a): {"reads": n, "breadth": b, "present": p, "virus": a}
            for s, a, n, b, p in rows}


def test_host_depletion_row():
    before = _long([("S1", "HERV", 100, 0.5, 1), ("S1", "PHAGE", 12, 0.02, 1),
                    ("S2", "HERV", 50, 0.4, 1)])
    after = _long([("S1", "PHAGE", 12, 0.02, 1)])
    row = v2lib.host_depletion_row(
        "S1", post_qc_pairs=1000,
        hostile={"reads_in": 2000, "reads_out": 1500, "reads_removed": 500},
        before=before, after=after, accession="HERV",
        pairs_before=60, pairs_after=6,
    )
    assert row["input_matches_post_qc"] == 1
    assert row["host_pct_of_post_qc"] == 25.0
    assert (row["herv_k113_reads_before"], row["herv_k113_reads_after"]) == (100, 0)
    # S2's HERV row must not leak into S1's count
    assert (row["refs_present_before"], row["refs_present_after"]) == (2, 1)


def test_host_depletion_flags_input_mismatch():
    row = v2lib.host_depletion_row(
        "S1", post_qc_pairs=1000,
        hostile={"reads_in": 1998, "reads_out": 1998, "reads_removed": 0},
        before={}, after={}, accession="X", pairs_before=0, pairs_after=0,
    )
    assert row["input_matches_post_qc"] == 0


# ---------------------------------------------------------------- agreement

def test_present_set_ignores_absent_rows():
    t = _long([("S1", "A", 20, 0.1, 1), ("S1", "B", 5, 0.5, 0), ("S2", "C", 30, 0.2, 1)])
    assert v2lib.present_set(t) == {"A", "C"}
    assert v2lib.present_set(t, "S1") == {"A"}


def test_agreement_partitions():
    ag = v2lib.agreement({"A", "B", "C"}, {"B", "C", "D"})
    assert ag == {"both": ["B", "C"], "only_a": ["A"], "only_b": ["D"]}


def test_compare_identical_tables_is_empty():
    t = _long([("S1", "A", 20, 0.1, 1)])
    assert v2lib.compare_long_tables(t, dict(t)) == []


def test_compare_catches_missing_row_and_breadth_change():
    old = _long([("S1", "A", 20, 0.1, 1), ("S1", "B", 11, 0.02, 1)])
    new = _long([("S1", "A", 20, 0.1001, 1)])
    keys = {(d["sample_id"], d["accession"]) for d in v2lib.compare_long_tables(old, new)}
    assert keys == {("S1", "A"), ("S1", "B")}


# ---------------------------------------------------------------- read classifier

def test_classify_unaligned_when_no_records():
    assert v2lib.classify_in_other("NC_1", [], 30) == "unaligned"


def test_classify_unaligned_when_only_unmapped_record():
    rec = sam(flag=0x1 | 0x4 | 0x40, rname="*", mapq=0)
    assert v2lib.classify_in_other("NC_1", [rec], 30) == "unaligned"


def test_classify_other_ref_uses_primary_not_secondary():
    # The secondary hit is on the target; the primary is elsewhere. The read
    # was placed elsewhere, so it must not be called same_ref.
    secondary = sam(flag=0x1 | 0x2 | 0x40 | 0x100, rname="NC_1", mapq=0)
    primary = sam(rname="NC_2", mapq=60)
    assert v2lib.classify_in_other("NC_1", [secondary, primary], 30) == "other_ref"


def test_classify_low_mapq_boundary():
    assert v2lib.classify_in_other("NC_1", [sam(mapq=29)], 30) == "same_ref_low_mapq"
    assert v2lib.classify_in_other("NC_1", [sam(mapq=30)], 30) == "same_ref_pass"


def test_classify_not_proper_pair():
    rec = sam(flag=0x1 | 0x40, mapq=60)
    assert v2lib.classify_in_other("NC_1", [rec], 30) == "same_ref_not_proper"


def test_passes_filters_matches_samtools_q_f2():
    assert v2lib.passes_filters(sam(mapq=30), 30)
    assert not v2lib.passes_filters(sam(mapq=29), 30)
    assert not v2lib.passes_filters(sam(flag=0x1 | 0x40), 30)
