"""Full per-sample profile of the leading candidate BioProject, plus paper hunt."""
import collections
import csv
import io
import re
import sys

sys.path.insert(0, "tools")
import ncbi  # noqa: E402

ACC = "PRJNA1183294"

_, csv_text = ncbi.runinfo(f"{ACC}[BioProject]")
rows = [r for r in csv.DictReader(io.StringIO(csv_text)) if r.get("Run")]
by_biosample = {r["BioSample"]: r for r in rows}
print(f"{len(rows)} runs / {len(by_biosample)} biosamples\n")

# Pull every BioSample record in batches.
attrs = {}
ids = sorted(by_biosample)
for i in range(0, len(ids), 20):
    xml = ncbi.efetch("biosample", ids[i : i + 20], "full", "xml")
    for block in re.findall(r"<BioSample .*?</BioSample>", xml, re.S):
        acc = re.search(r'accession="([^"]+)"', block)
        if not acc:
            continue
        d = {}
        for k, v in re.findall(
            r'<Attribute attribute_name="([^"]+)"[^>]*>(.*?)</Attribute>', block, re.S
        ):
            d[k] = re.sub(r"\s+", " ", v).strip()
        attrs[acc.group(1)] = d

print("attribute keys present across samples:")
keycount = collections.Counter(k for d in attrs.values() for k in d)
for k, n in keycount.most_common():
    print(f"   {k}: {n}/{len(attrs)}")

print("\nDiagnosis distribution:")
for val, n in collections.Counter(d.get("Diagnosis", "<absent>") for d in attrs.values()).most_common():
    print(f"   {n:3d}  {val}")

print("\ngeo_loc_name distribution:")
for val, n in collections.Counter(d.get("geo_loc_name", "<absent>") for d in attrs.values()).most_common():
    print(f"   {n:3d}  {val}")

print("\nisolation_source distribution:")
for val, n in collections.Counter(d.get("isolation_source", "<absent>") for d in attrs.values()).most_common():
    print(f"   {n:3d}  {val}")

print("\nDiagnosis x country:")
cross = collections.Counter(
    (d.get("geo_loc_name", "?"), d.get("Diagnosis", "?")) for d in attrs.values()
)
for (c, dg), n in sorted(cross.items()):
    print(f"   {c:<12} {dg:<28} {n}")

print("\nper-sample table (biosample, run, name, diagnosis, country, spots, MB):")
for bsa in ids:
    r = by_biosample[bsa]
    d = attrs.get(bsa, {})
    print(
        f"   {bsa}  {r['Run']}  {r['SampleName']:<10} "
        f"{d.get('Diagnosis', '?'):<26} {d.get('geo_loc_name', '?'):<10} "
        f"{int(r['spots']):>11,}  {r['size_MB']:>6} MB"
    )

print("\n\n########## paper search")
for term in [
    "DENTAID AND subgingival AND metagenom*",
    "subgingival microbiome AND periodontitis AND four countries",
    "(periodontitis[Title]) AND (metagenom*[Title]) AND (countries[Title/Abstract])",
    "PRJNA1183294",
]:
    count, pmids = ncbi.esearch("pubmed", term, retmax=6)
    print(f"\n-- {term!r}: {count} hits")
    for rec in ncbi.esummary("pubmed", pmids[:6]):
        print(f"   PMID {rec.get('uid')}: {rec.get('title', '')[:120]}")
        print(f"        {rec.get('fulljournalname', '')} {rec.get('pubdate', '')}")
