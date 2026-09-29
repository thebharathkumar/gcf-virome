# Metadata standards used in this project

Every format rule the validator enforces is recorded here with the source it
came from. Nothing in this file is from memory. Where a check is looser than
the real standard, that is stated.

Verified on 2026-09-27.

## 1. REDCap data dictionary CSV

The data dictionary has 18 columns in a fixed order. The header row below was
taken from a real REDCap-exported data dictionary, not from documentation
prose:

```
"Variable / Field Name","Form Name","Section Header","Field Type","Field Label","Choices, Calculations, OR Slider Labels","Field Note","Text Validation Type OR Show Slider Number","Text Validation Min","Text Validation Max",Identifier?,"Branching Logic (Show field only if...)","Required Field?","Custom Alignment","Question Number (surveys only)","Matrix Group Name","Matrix Ranking?","Field Annotation"
```

Points that are easy to get wrong and are therefore checked:

- `Identifier?`, `Required Field?` and `Matrix Ranking?` all carry a question
  mark. Dropping it is a common error and the validator rejects it.
- `Branching Logic (Show field only if...)` contains three literal periods.
- Mandatory columns are A, B, D and E (`Variable / Field Name`, `Form Name`,
  `Field Type`, `Field Label`). Column F is mandatory only for
  dropdown, radio, checkbox and calc fields.

| Setting | Rule enforced |
|---|---|
| Variable names | lowercase letters, digits, underscores; start with a letter; 2 to 26 characters |
| Field types | text, notes, calc, dropdown, radio, checkbox, yesno, truefalse, file, slider, descriptive, sql |
| Choice syntax | `1, Label One \| 2, Label Two`, code then comma then label, choices separated by a pipe |
| Validation types | lowercase, for example `date_ymd`, not `DATE_YMD` |
| Required marker | lowercase `y`, or empty. `Y`, `TRUE` and `yes` are rejected |

Sources:
- <https://raw.githubusercontent.com/redcap-tools/redcap-test-datasets/master/case-01/test-case-01-data-dictionary.csv>
- <https://cri.uchicago.edu/wp-content/uploads/2015/12/REDCap-Data-Dictionary.pdf>
- <https://www.iths.org/wp-content/uploads/201-Data-Dictionary.pdf>

Not verified: whether REDCap enforces a hard maximum variable-name length.
Sources disagree between 26 recommended and 100 maximum. The validator uses
26, which is the stricter reading. The `sql` field type is accepted but its
admin-only status was not confirmed from a primary source.

## 2. REDCap records CSV

- The first column must be the record ID field, which is the first field of
  the first form in the dictionary.
- Column headers must match dictionary variable names exactly.
- Checkbox fields appear as one column per choice, named
  `variable___code` with three underscores, each holding `0` or `1`.
- Radio, dropdown, yesno and truefalse fields hold the raw numeric code, not
  the display label. Writing `Belgium` where `1` is expected is one of the
  most common import failures and is tested for.
- `date_ymd` fields are written `YYYY-MM-DD`.

Sources:
- <https://www.iths.org/wp-content/uploads/REDCap-Importing-Exporting-302.pdf>
- <https://guides.temple.edu/c.php?g=936400&p=8404930>

## 3. NCBI BioSample package

Package used: **MIMS.me.human-oral.6.0**, "MIMS: metagenome/environmental,
human-oral; version 6.0". This is the most specific package for human oral
metagenome samples. The generic `MIMS.me.host-associated.6.0` has the same
mandatory set and would also be acceptable.

Seven mandatory package attributes, plus the `sample_name` and `organism`
base columns:

```
collection_date  env_broad_scale  env_local_scale  env_medium
geo_loc_name     host             lat_lon
```

Format rules enforced:

| Attribute | Rule |
|---|---|
| `collection_date` | `YYYY`, `YYYY-MM`, `YYYY-MM-DD` or ISO 8601. A range is two dates separated by `/` |
| `geo_loc_name` | INSDC country name, optionally `Country:Region` |
| `lat_lon` | decimal degrees with compass letters, `50.8802 N 4.6934 E`. Signed decimals such as `50.88, 4.69` are rejected |
| `env_*` | term text with an ontology ID in brackets |
| missing values | `not collected`, `not applicable` or `missing` |

Sources:
- <https://www.ncbi.nlm.nih.gov/biosample/docs/packages/>
- <https://www.ncbi.nlm.nih.gov/biosample/docs/packages/MIMS.me.human-oral.6.0/>
- <https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/>
- <https://www.ncbi.nlm.nih.gov/biosample/docs/submission/faq/>

Not verified: the literal template header row with leading asterisks on
mandatory fields. NCBI's template generator is a JavaScript application and
the raw file could not be downloaded. The attribute names themselves come
from the package page and are reliable.

## 4. The environmental triplet

This is the part of the submission that is most often filled in carelessly,
so each term was resolved against the EBI Ontology Lookup Service and checked
for obsolescence before use.

Every ID below links to its Ontology Lookup Service entry. Click it and the
`is_obsolete` flag and label are on the page.

| Field | Term used | Status | Verify |
|---|---|---|---|
| `env_broad_scale` | `oral cavity [UBERON:0000167]` | exists, not obsolete | [OLS](https://www.ebi.ac.uk/ols4/ontologies/uberon/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FUBERON_0000167) |
| `env_local_scale` | `gingiva [UBERON:0001828]` | exists, not obsolete | [OLS](https://www.ebi.ac.uk/ols4/ontologies/uberon/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FUBERON_0001828) |
| `env_medium` | `bodily fluid material [ENVO:02000019]` | exists, not obsolete | [OLS](https://www.ebi.ac.uk/ols4/ontologies/envo/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FENVO_02000019) |

Two plausible-looking choices were rejected after checking:

| Rejected ID | What it actually is | Verify |
|---|---|---|
| `ENVO:00009003` | **"obsolete human-associated habitat"**, `is_obsolete: true`, no recorded replacement. Widely used in real submissions regardless. | [OLS](https://www.ebi.ac.uk/ols4/ontologies/envo/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FENVO_00009003) |
| `UBERON:0006932` | **"vestibular epithelium"**, not gingival crevicular fluid. | [OLS](https://www.ebi.ac.uk/ols4/ontologies/uberon/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FUBERON_0006932) |

To re-check any of these from the command line, which is the exact call used
when the terms were chosen:

```bash
curl -sS -G "https://www.ebi.ac.uk/ols4/api/ontologies/envo/terms" \
  --data-urlencode "iri=http://purl.obolibrary.org/obo/ENVO_00009003" \
  | python3 -c "import sys,json; d=json.load(sys.stdin)['_embedded']['terms'][0]; \
print(d['label'], '| is_obsolete:', d['is_obsolete'])"
# obsolete human-associated habitat | is_obsolete: True
```

There is **no ENVO or UBERON term for gingival crevicular fluid**. The nearest
matches in any ontology are `BTO:0003364` "gingival fluid" and MeSH
`D005883`. Because MIxS expects ENVO for `env_medium`, the nearest correct
ENVO parent is used and the precise specimen is carried in
`isolation_source` and `samp_collect_device` instead. This is a deliberate
curation decision, not an oversight.

MIxS recommends that `env_broad_scale` be a subclass of `biome [ENVO:00000428]`
and `env_medium` a subclass of `environmental material [ENVO:00010483]`. Both
parent terms were confirmed to exist. No ENVO biome subclass describes a human
body site, which is why a UBERON anatomical term is used for the broad scale.
A curator submitting this for real should expect NCBI to accept it but may
wish to confirm with the BioSample help desk.

## 5. SRA library and run sheet

NCBI keeps sample attributes and sequencing metadata in two separate sheets.
Column names below are from NCBI's own SRA metadata spreadsheet
(`SRA_metadata_acc_example.xlsx`), parsed directly:

```
bioproject_accession  biosample_accession  library_ID  title
library_strategy  library_source  library_selection  library_layout
platform  instrument_model  design_description  filetype
filename  filename2  filename3  filename4
```

Controlled vocabularies enforced, taken from the "Library and Platform Terms"
sheet of that file:

- `library_strategy`: 24 terms. This dataset uses `WGS`.
- `library_source`: 7 terms. This dataset uses `METAGENOMIC`.
- `library_selection`: 27 terms. This dataset uses `size fractionation`, which
  is the value the original submitters recorded in SRA.
- `library_layout`: `single` or `paired`, lowercase.
- `platform`: 9 terms. This dataset uses `ILLUMINA`.
- `filetype`: 9 terms. This dataset uses `fastq`.

Not verified: `Illumina NovaSeq 6000` does not appear in the instrument list
of the example spreadsheet NCBI links from its documentation, whose Illumina
list stops around the HiSeq 4000 and NextSeq 550 era. The file appears to be
stale. NovaSeq 6000 is one of the most common instrument values in current SRA
records and is the value in this BioProject's own SRA entries, so it is used
here, but it was not confirmed against a current template. The validator does
not check `instrument_model` against a fixed list for this reason.

`filename` and `filename2` name the subsampled FASTQ files this pipeline
produces. They are illustrative. Nothing here is submitted to NCBI.
