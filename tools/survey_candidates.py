"""List BioProjects matching oral/salivary virome queries, with descriptions."""
import sys

sys.path.insert(0, "tools")
import ncbi  # noqa: E402

QUERIES = [
    (
        "oral/saliva + viral terms",
        "(oral[Title] OR saliva[Title] OR salivary[Title] OR subgingival[Title] "
        "OR periodontitis[Title] OR periodontal[Title] OR dental[Title]) "
        "AND (virome[Title] OR virus[Title] OR viral[Title] OR phage[Title])",
    ),
    (
        "periodontitis + metagenome",
        "(periodontitis[Title] OR periodontal[Title] OR caries[Title]) "
        "AND (metagenom*[Title] OR microbiome[Title] OR shotgun[Title])",
    ),
]

seen = {}
for label, term in QUERIES:
    count, ids = ncbi.esearch("bioproject", term, retmax=60)
    print(f"\n########## {label}: {count} hits, showing {len(ids)}\n")
    for rec in ncbi.esummary("bioproject", ids):
        acc = rec.get("project_acc", "?")
        if acc in seen:
            continue
        seen[acc] = True
        title = (rec.get("project_title") or "").strip()
        desc = (rec.get("project_description") or "").strip().replace("\n", " ")
        org = (rec.get("organism_name") or "").strip()
        print(f"{acc}  [{org}]  {title[:130]}")
        if desc:
            print(f"        {desc[:260]}")
