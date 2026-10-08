"""Pure functions shared by the results/v2 experiments.

Kept free of Snakemake and subprocess calls so they can be unit-tested on
hand-built inputs. The scripts that Snakemake runs do the file and process
handling and call into here for every number they report.
"""
import csv
import json
import re

# A read counts toward a reference only if it passes the same filters as the
# original run: MAPQ >= min_mapq and the SAM proper-pair flag.
FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_UNMAP = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def flagstat_proper_pairs(text):
    """Properly paired reads from samtools flagstat, halved to pairs.

    Identical to the rule used by read_retention.py in the original run, so a
    "viral pairs" number here means the same thing as mapped_pairs there.
    """
    m = re.search(r"^(\d+) \+ \d+ properly paired", text, re.M)
    if m:
        return int(m.group(1)) // 2
    m = re.search(r"^(\d+) \+ \d+ in total", text, re.M)
    return (int(m.group(1)) // 2) if m else 0


def fastp_post_qc_pairs(fastp_json_text):
    d = json.loads(fastp_json_text)
    return d["summary"]["after_filtering"]["total_reads"] // 2


def hostile_stats(hostile_json_text):
    """Hostile writes a JSON list with one dict per input; we pass one pair."""
    d = json.loads(hostile_json_text)
    if isinstance(d, list):
        if len(d) != 1:
            raise ValueError(f"expected one Hostile record, found {len(d)}")
        d = d[0]
    return {
        "reads_in": int(d["reads_in"]),
        "reads_out": int(d["reads_out"]),
        "reads_removed": int(d["reads_removed"]),
    }


def read_long_table(path):
    """virus_counts_long.tsv -> {(sample, accession): row dict}."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["sample_id"], row["accession"])] = {
                "reads": int(row["reads"]),
                "breadth": float(row["breadth"]),
                "present": int(row["present"]),
                "virus": row["virus"],
            }
    return out


# ---------------------------------------------------------------------------
# Host depletion
# ---------------------------------------------------------------------------

def host_depletion_row(sample, post_qc_pairs, hostile, before, after,
                       accession, pairs_before, pairs_after):
    """One per-sample row of the host depletion comparison.

    before / after are long tables from read_long_table for the same sample
    set. accession is the HERV-K113 RefSeq accession.
    """
    reads_in = hostile["reads_in"]
    removed = hostile["reads_removed"]
    herv_before = before.get((sample, accession), {}).get("reads", 0)
    herv_after = after.get((sample, accession), {}).get("reads", 0)
    present_before = sum(
        1 for (s, _), r in before.items() if s == sample and r["present"]
    )
    present_after = sum(
        1 for (s, _), r in after.items() if s == sample and r["present"]
    )
    return {
        "sample_id": sample,
        "post_qc_pairs": post_qc_pairs,
        "hostile_reads_in": reads_in,
        # Hostile counts reads, fastp pairs are reads / 2. If these disagree
        # the depletion step did not see the same input as the original run.
        "input_matches_post_qc": int(reads_in == 2 * post_qc_pairs),
        "host_reads_removed": removed,
        "host_pct_of_post_qc": (
            round(100.0 * removed / (2 * post_qc_pairs), 4) if post_qc_pairs else 0.0
        ),
        "viral_pairs_before": pairs_before,
        "viral_pairs_after": pairs_after,
        "herv_k113_reads_before": herv_before,
        "herv_k113_reads_after": herv_after,
        "refs_present_before": present_before,
        "refs_present_after": present_after,
    }


# ---------------------------------------------------------------------------
# Aligner agreement
# ---------------------------------------------------------------------------

def present_set(long_table, sample=None):
    """Accessions passing the presence filter, optionally within one sample."""
    return {
        acc for (s, acc), r in long_table.items()
        if r["present"] and (sample is None or s == sample)
    }


def agreement(a, b):
    a, b = set(a), set(b)
    return {"both": sorted(a & b), "only_a": sorted(a - b), "only_b": sorted(b - a)}


def compare_long_tables(old, new):
    """Rows where two runs of the same step disagree. Empty means identical."""
    diffs = []
    for key in sorted(set(old) | set(new)):
        o, n = old.get(key), new.get(key)
        o_reads = o["reads"] if o else 0
        n_reads = n["reads"] if n else 0
        o_pres = o["present"] if o else 0
        n_pres = n["present"] if n else 0
        o_br = o["breadth"] if o else 0.0
        n_br = n["breadth"] if n else 0.0
        if (o_reads, o_pres) != (n_reads, n_pres) or abs(o_br - n_br) > 1e-6:
            diffs.append({
                "sample_id": key[0], "accession": key[1],
                "reads_v1": o_reads, "reads_v2": n_reads,
                "breadth_v1": o_br, "breadth_v2": n_br,
                "present_v1": o_pres, "present_v2": n_pres,
            })
    return diffs


# ---------------------------------------------------------------------------
# Read-level explanation of aligner disagreement
# ---------------------------------------------------------------------------

def parse_sam_line(line):
    f = line.rstrip("\n").split("\t")
    tags = {}
    for t in f[11:]:
        k, _, v = t.split(":", 2)
        tags[k] = v
    cigar = f[5]
    soft = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op == "S")
    return {
        "qname": f[0],
        "flag": int(f[1]),
        "rname": f[2],
        "pos": int(f[3]),
        "mapq": int(f[4]),
        "cigar": cigar,
        "soft_clipped": soft,
        "nm": int(tags["NM"]) if "NM" in tags else None,
        "mate": 2 if int(f[1]) & 0x80 else 1,
    }


def passes_filters(rec, min_mapq):
    return (
        not rec["flag"] & FLAG_UNMAP
        and rec["flag"] & FLAG_PROPER
        and rec["mapq"] >= min_mapq
    )


def classify_in_other(target_ref, other_records, min_mapq):
    """Why a read counted on target_ref by one aligner was not by the other.

    other_records: parsed SAM records the other aligner wrote for the same
    read and mate (primary, secondary and supplementary). Returns one label.

      unaligned           the other aligner reported no alignment at all
      same_ref_pass       same reference and passes filters (difference is
                          in read count or breadth, not in this read)
      same_ref_low_mapq   same reference, MAPQ below the cutoff
      same_ref_not_proper same reference, MAPQ ok, not a proper pair
      other_ref           primary alignment is on a different reference
    """
    mapped = [r for r in other_records if not r["flag"] & FLAG_UNMAP]
    if not mapped:
        return "unaligned"
    primary = [r for r in mapped if not r["flag"] & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)]
    rec = primary[0] if primary else mapped[0]
    if rec["rname"] != target_ref:
        return "other_ref"
    if rec["mapq"] < min_mapq:
        return "same_ref_low_mapq"
    if not rec["flag"] & FLAG_PROPER:
        return "same_ref_not_proper"
    return "same_ref_pass"
