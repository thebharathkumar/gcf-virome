import sys

sys.path.insert(0, "tools")
import ncbi  # noqa: E402

for pmid in ["40202358", "41522492"]:
    print("=" * 100)
    print(ncbi.efetch("pubmed", [pmid], "abstract", "text")[:4000])
