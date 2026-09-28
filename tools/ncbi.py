"""Thin NCBI E-utilities helper used to verify accessions before any download.

Read-only. Every accession quoted in the write-up must survive a call here.
"""
import json
import sys
import time
import urllib.parse
import urllib.request

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(endpoint, params):
    url = f"{E}/{endpoint}?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def esearch(db, term, retmax=100):
    raw = _get("esearch.fcgi", {"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    res = json.loads(raw)["esearchresult"]
    return int(res["count"]), res.get("idlist", [])


def esummary(db, ids):
    if not ids:
        return []
    raw = _get("esummary.fcgi", {"db": db, "id": ",".join(ids), "retmode": "json"})
    res = json.loads(raw)["result"]
    return [res[u] for u in res["uids"]]


def efetch(db, ids, rettype, retmode="text"):
    return _get(
        "efetch.fcgi",
        {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": retmode},
    )


# NCBI sometimes returns the runinfo CSV with no header line, so carry the
# canonical column order and re-attach it when it is absent.
RUNINFO_HEADER = (
    "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,"
    "AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,"
    "LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,"
    "Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,"
    "SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,"
    "g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,"
    "Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,"
    "dbgap_study_accession,Consent,RunHash,ReadHash"
)


def runinfo(term):
    """SRA run table for a search term, as CSV text with a guaranteed header."""
    count, ids = esearch("sra", term, retmax=1000)
    if not ids:
        return count, ""
    text = efetch("sra", ids, "runinfo", "csv").lstrip("\n")
    if not text.startswith("Run,"):
        text = RUNINFO_HEADER + "\n" + text
    return count, text


def elink(dbfrom, db, ids):
    raw = _get(
        "elink.fcgi",
        {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "json"},
    )
    out = []
    for link_set in json.loads(raw).get("linksets", []):
        for ldb in link_set.get("linksetdbs", []):
            out.extend(ldb.get("links", []))
    return sorted(set(out))
