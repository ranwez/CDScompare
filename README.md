# CDScompare

**CDScompare performs the quantitative comparison of structural genome annotations (GFF3) of the same genome.**  
For each pair of overlapping genes, it computes numerical similarity scores based on CDS organization, exon–intron structure and reading frame conservation.  
It is implemented in Python and distributed as a command-line tool.


The tool supports:
- pairwise comparison between two genome annotations
- multi-comparison of several annotations against a common reference
- two pairing strategies (best or all)

## Installation

### Requirements
- Python ≥ 3.9 with `pip` 

### Install from the Git repository

```bash
git clone git@github.com:ranwez/CDScompare.git
cd CDScompare
pip install .
```
*NB: Once installed, the repository is no longer required.*  

Successful installation can be tested with:
```bash
cdscompare --version
```


## Command-line interface

```bash
cdscompare [OPTIONS] GFF_1 GFF_2 [GFF_3 ...]
```

### Positional arguments

| Argument | Description |
|---------|-------------|
| `GFF_1` | Annotation file in GFF3 format. |
| `GFF_2 ...` | One or more additional annotation files. When more than two GFF files are provided, the first file (`GFF_1`) is used as the reference, and all other files are compared against it. |

> **Note:** 
> - At least **two GFF files** are required.
> - In pairwise comparisons (exactly two GFF files), identity scores and gene pairings are invariant to the input file order.
> - GFF input file basenames must be unique (used as annotation identifiers).


### Options

| Option | Description |
|-------|-------------|
| `-d, --out_dir` | Output directory where result files are written (default: `results`). |
| `-p, --pairing_mode` | Pairing strategy used within clusters of overlapping genes. Possible values are:<br>• `best` (default): selects a globally optimal gene pairing using dynamic programming.<br>• `all`: reports all overlapping gene pairings without global optimization. |


## Output files

### Pairwise comparison

CDScompare will produce:

#### 1. Detailed comparison file (`<annotation1>_vs_<annotation2>.csv`)

This file contains one line per gene comparison, with the following columns:

| Column | Description |
|------|-------------|
| `chromosome` | Chromosome identifier including strand (`_direct` or `_reverse`). |
| `cluster` | Identifier of the cluster of overlapping genes. |
| `annot1_gene` | Gene identifier in the first input annotation. |
| `annot2_gene` | Gene identifier in the second input annotation. |
| `matches` | Number of nucleotide positions matching between CDS structures. |
| `mismatches` | Total number of mismatched nucleotides (exon–intron + reading frame mismatches). |
| `identity_score` | Percentage identity computed as matches / (matches + mismatches). |
| `annot1_start` | Genomic start coordinate of the reference gene. |
| `annot1_end` | Genomic end coordinate of the reference gene. |
| `annot2_start` | Genomic start coordinate of the alternative gene. |
| `annot2_end` | Genomic end coordinate of the alternative gene. |
| `annot1_mRNA` | Identifier of the selected reference mRNA used for the comparison. |
| `annot2_mRNA` | Identifier of the selected alternative mRNA used for the comparison. |
| `EI_mismatches_zones` | Genomic intervals where exon–intron structures differ between annotations. |
| `RF_mismatches_zones` | Genomic intervals where reading frames differ between annotations. |
| `EI_mismatches` | Total length (in nucleotides) of exon–intron mismatches. |
| `RF_mismatches` | Total length (in nucleotides) of reading frame mismatches. |
| `annot1_mRNA_number` | Number of mRNAs annotated for the reference gene. |
| `annot2_mRNA_number` | Number of mRNAs annotated for the alternative gene. |

Special values:
- `_` : undefined or not applicable
- `~` : no corresponding gene reported for this comparison. This occurs when:
  - no overlapping gene exists in the other annotation, or
  - the gene was paired with a different gene in `best` pairing mode.


#### 2. Summary file (`<annotation1>_vs_<annotation2>.txt`)

This file contains a global summary:

```
Number of loci (whole data):
- found in both annotations : X
- found only in the reference : Y
- found only in the alternative : Z
```

### Multi-comparison mode

When more than two annotation files are provided:

- A detailed comparison file and a summary file are produced for each alternative annotation.

- An additional global summary file (`synthesis_<annotation1>.csv`) is generated, with one line per reference annotation gene, and two columns for each alternative annotation:

| Column | Description |
|--------|---------|
| `Reference_locus` | Gene identifier in the reference annotation. |
| `<alt>.locus` | Identifier of the best matching gene in the corresponding alternative annotation. |
| `<alt>.identity` | Identity score (%) for the corresponding gene pairing. |


## Scope and assumptions

- Gene overlap is detected based on **gene genomic coordinates**.

- Identity scores are computed **only from CDS features**.
  UTRs are ignored, and exon–intron structure is inferred from CDS organization.

- Alternative splicing is handled by selecting the **best matching mRNA pair** for each gene pairing.

- CDScompare assumes that all input annotation files:
  - describe the **same genome assembly**
  - follow standard **GFF3 conventions**, with the following feature hierarchy:
    - `gene` features with an `ID` attribute
    - `mRNA` features with `ID` and `Parent` attributes
    - `CDS` features with a `Parent` attribute pointing to an mRNA
  - contain **non-overlapping CDS coordinates within a single mRNA**
  - use **unique gene identifiers (`ID`) within each file**



## Citation

...


## Developers notes

### Development installation

This method installs the package in editable mode and includes development dependencies (pytest and ruff).
```bash
git clone git@github.com:ranwez/CDScompare.git
cd CDScompare
pip install -e .[dev]
```

To run the tests:
```bash
python -m pytest -v
```

### Developer documentation
The internal code structure documentation is generated using Sphinx.

To build the documentation locally:
```bash
pip install -e .[docs]
python -m sphinx -b html docs/source docs/build/html
```
The generated documentation can be opened locally in a web browser (docs/build/html/index.html).
