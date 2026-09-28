"""Run build_metadata.py and the validator outside Snakemake.

The Snakemake script directive injects a `snakemake` global. This harness
fakes it so the metadata scripts can be exercised while the main workflow is
still downloading, rather than discovering a bug two hours from now.
"""
import runpy
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build" / "smoke_metadata"
OUT.mkdir(parents=True, exist_ok=True)

import yaml  # noqa: E402

cfg = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())

fake = types.SimpleNamespace(
    input=types.SimpleNamespace(samples=str(ROOT / "config" / "samples.tsv")),
    output=types.SimpleNamespace(
        dictionary=str(OUT / "redcap_data_dictionary.csv"),
        records=str(OUT / "redcap_records.csv"),
        biosample=str(OUT / "biosample_MIMS_human_oral.tsv"),
        sra=str(OUT / "sra_metadata.tsv"),
    ),
    config=cfg,
)

import builtins  # noqa: E402

builtins.snakemake = fake
runpy.run_path(str(ROOT / "workflow" / "scripts" / "build_metadata.py"), run_name="__main__")

print("\n--- validating ---\n")
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
import validate_metadata as vm  # noqa: E402

code = vm.main([
    "--dictionary", fake.output.dictionary,
    "--records", fake.output.records,
    "--biosample", fake.output.biosample,
    "--sra", fake.output.sra,
    "--manifest", fake.input.samples,
    "--report", str(OUT / "validation_report.txt"),
])
print(f"exit code: {code}")
sys.exit(code)
