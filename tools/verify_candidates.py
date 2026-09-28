"""Verify shortlisted BioProjects: run counts, strategy, data volume, metadata.

Nothing is downloaded here. This only reads NCBI metadata so that every
accession quoted in the write-up is known to exist and known to be usable.
"""
import collections
import csv
import io
import re
import sys

sys.path.insert(0, "tools")
import ncbi  # noqa: E402

CANDIDATES = [
    ("PRJNA1183294", "Subgingival shotgun metagenome, healthy vs periodontitis, 4 countries"),
    ("PRJNA663926", "Saliva virome+microbiome, paediatric HFMD vs healthy controls"),
    ("PRJNA1223632", "Whole salivary virome, Sjogren disease vs healthy controls"),
    ("PRJNA692713", "Oral microbiome+virome, OSCC vs normal, shotgun"),
    ("PRJNA1265940", "Oral cavity WGS metagenome, mild vs severe periodontitis"),
    ("PRJNA1238390", "Subgingival metagenome, health vs grades of periodontitis"),
    ("PRJNA1482867", "Dental plaque shotgun, before/after periodontal therapy"),
]

SKIP_VALUES = {"not applicable", "missing", "not collected", "na", "n/a", "none", ""}


def summarise(acc):
    _, csv_text = ncbi.runinfo(f"{acc}[BioProject]")
    if not csv_text.strip():
        print("  no SRA runs found")
        return []
    rows = [r for r in csv.DictReader(io.StringIO(csv_text)) if r.get("Run")]
    if not rows:
        print("  runinfo returned no parsable rows")
        return []

    def counts(field):
        return dict(collections.Counter(r[field] for r in rows if r.get(field)))

    spots = [int(r["spots"]) for r in rows if r.get("spots", "").isdigit()]
    bases = [int(r["bases"]) for r in rows if r.get("bases", "").isdigit()]
    mb = [int(r["size_MB"]) for r in rows if r.get("size_MB", "").isdigit()]

    print(f"  runs={len(rows)}  biosamples={len({r['BioSample'] for r in rows})}")
    print(f"  strategy={counts('LibraryStrategy')}  selection={counts('LibrarySelection')}")
    print(f"  source={counts('LibrarySource')}  layout={counts('LibraryLayout')}")
    print(f"  model={counts('Model')}")
    print(f"  organism={counts('ScientificName')}")
    if spots:
        print(f"  spots/run: median={sorted(spots)[len(spots) // 2]:,}  total={sum(spots):,}")
    if bases:
        print(f"  bases: total={sum(bases) / 1e9:.1f} Gbp  median/run={sorted(bases)[len(bases) // 2] / 1e6:,.0f} Mbp")
    if mb:
        print(f"  size_MB: total={sum(mb):,} MB  median/run={sorted(mb)[len(mb) // 2]:,} MB")
    print(f"  example runs: {[r['Run'] for r in rows[:3]]}")
    print(f"  sample names: {[r['SampleName'] for r in rows[:6]]}")
    return sorted({r["BioSample"] for r in rows})


def biosample_attributes(biosamples, n=4):
    if not biosamples:
        return
    xml = ncbi.efetch("biosample", biosamples[:n], "full", "xml")
    for block in re.findall(r"<BioSample .*?</BioSample>", xml, re.S)[:n]:
        acc = re.search(r'accession="([^"]+)"', block)
        title = re.search(r"<Title>(.*?)</Title>", block, re.S)
        print(f"    --- {acc.group(1) if acc else '?'} :: {(title.group(1) if title else '').strip()[:70]}")
        for k, v in re.findall(
            r'<Attribute attribute_name="([^"]+)"[^>]*>(.*?)</Attribute>', block, re.S
        ):
            v = re.sub(r"\s+", " ", v).strip()
            if v.lower() in SKIP_VALUES:
                continue
            print(f"        {k} = {v[:90]}")


def linked_papers(acc):
    _, ids = ncbi.esearch("bioproject", f"{acc}[Project Accession]", retmax=5)
    if not ids:
        return
    pmids = ncbi.elink("bioproject", "pubmed", ids)
    if not pmids:
        print("    paper: no PubMed link recorded in BioProject")
        return
    for rec in ncbi.esummary("pubmed", pmids[:4]):
        print(f"    paper: PMID {rec.get('uid')} :: {rec.get('title', '')[:110]}")
        print(f"            {rec.get('fulljournalname', '')} {rec.get('pubdate', '')}")


for acc, note in CANDIDATES:
    print(f"\n================ {acc}  {note}")
    try:
        bs = summarise(acc)
        linked_papers(acc)
        biosample_attributes(bs)
    except Exception as exc:  # noqa: BLE001
        print(f"  ERROR: {type(exc).__name__}: {exc}")
