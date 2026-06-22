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

One of the following must be available:

- Python ≥ 3.9 with `pip` or `pipx`
- Docker or Apptainer / Singularity

---

### Installation from PyPI (recommended)

#### Using `pipx` (isolated CLI install)

```bash
pipx install cdscompare==<version>
```

_Replace \<version\> with a specific release number (e.g. 0.3.0rc5)._

#### Using `pip`

```bash
pip install cdscompare==<version>
```

Test the installation:

```bash
cdscompare --version
```

### Container image (Docker / Apptainer)

A pre-built OCI image is available via GitHub Container Registry:

```bash
ghcr.io/johgi/cdscompare:<version>
```

#### Docker

Pull the image:

```bash
docker pull ghcr.io/johgi/cdscompare:<version>
```

To run on local files:

```bash
docker run --rm \
  -v "$PWD:/work" \
  -w /work \
  ghcr.io/johgi/cdscompare:<version> \
  file1.gff file2.gff
```

#### Apptainer / Singularity

Pull the image:

```bash
apptainer pull docker://ghcr.io/johgi/cdscompare:<version>
```

To run on local files:

```bash
apptainer run \
  --bind "$PWD:/work" \
  cdscompare_<version>.sif \
  /work/file1.gff /work/file2.gff
```

### Installing development / unreleased versions

```bash
git clone git@github.com:ranwez/CDScompare.git
cd CDScompare
pip install .
```

_After installation, the repository is no longer required._

---

## Preparing input GFF files

CDScompare expects both input annotation files to describe protein-coding gene models in a consistent GFF3-like structure. Because GFF/GTF files may differ substantially between annotation tools and databases, we strongly recommend standardizing input files before running CDScompare.

In particular, each input file should contain:

- tab-delimited GFF3 records with 9 columns;
- `gene`, `mRNA` and `CDS` features organized in a consistent `gene -> mRNA -> CDS` hierarchy;
- valid and unique `ID`/`Parent` attributes linking genes, transcripts and CDS features;
- coherent genomic coordinates;
- non-overlapping CDS features within each transcript;
- valid CDS phases (`0`, `1`, `2`, or `.`);
- compatible seqids in column 1 between the two annotations.

The two annotations must refer to the same genome assembly and coordinate system.

CDScompare compares only coding sequences. Non-coding genes, pseudogenes without CDS features, and transcripts without CDS features are not considered.

### Seqid compatibility

The seqids in column 1 of the two GFF files must be compatible. CDScompare can harmonize an optional leading `chr` prefix when this improves the number of shared seqids and does not merge distinct seqids within the same annotation. For example, annotations using `chr1`, `chr2`, `chr3` can therefore be compared to annotations using `1`, `2`, `3`.

More complex naming differences, such as `NC_000001.11` vs `chr1`, `chromosome_1` vs `1`, or scaffold-specific naming schemes, must be harmonized by the user before running CDScompare.

### Recommended preprocessing with AGAT

Before running CDScompare, we strongly recommend standardizing each input annotation with [AGAT](https://github.com/NBISweden/AGAT) (Dainat, 2022). For example:

```bash
agat_convert_sp_gxf2gxf.pl -g annotation.gff3 -o annotation.agat.gff3
```

AGAT can help standardize GFF/GTF files, sort features, fix duplicated IDs, add missing `ID` and `Parent` attributes, add missing parent features when possible, and produce a more consistent GFF3 structure. This preprocessing step is especially recommended when input files come from different annotation tools or databases.

---

## Command-line interface

```bash
cdscompare [OPTIONS] GFF_1 GFF_2 [GFF_3 ...]
```

### Positional arguments

| Argument    | Description                                                                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `GFF_1`     | Annotation file in GFF3 format.                                                                                                                                                     |
| `GFF_2 ...` | One or more additional annotation files. When more than two GFF files are provided, the first file (`GFF_1`) is used as the reference, and all other files are compared against it. |

> **Note:**
>
> - At least **two GFF files** are required.
> - In pairwise comparisons (exactly two GFF files), identity scores and gene pairings are invariant to the input file order.
> - GFF input file basenames must be unique (used as annotation identifiers).

### Options

| Option               | Description                                                                                                                                                                                                                                               |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-d, --out_dir`      | Output directory where result files are written (default: `results`).                                                                                                                                                                                     |
| `-p, --pairing_mode` | Pairing strategy used within clusters of overlapping genes. Possible values are:<br>• `best` (default): selects a globally optimal gene pairing using dynamic programming.<br>• `all`: reports all overlapping gene pairings without global optimization. |

## Output files

### Pairwise comparison

CDScompare will produce:

#### 1. Detailed comparison file (`<annotation1>_vs_<annotation2>.csv`)

This file contains one line per gene comparison, with the following columns:

| Column                | Description                                                                      |
| --------------------- | -------------------------------------------------------------------------------- |
| `chromosome`          | Chromosome identifier including strand (`_direct` or `_reverse`).                |
| `cluster`             | Identifier of the cluster of overlapping genes.                                  |
| `annot1_gene`         | Gene identifier in the first input annotation.                                   |
| `annot2_gene`         | Gene identifier in the second input annotation.                                  |
| `matches`             | Number of nucleotide positions that are coding in both annotations and in the same reading frame.                  |
| `mismatches`          | Total number of mismatched nucleotides (coding/non-coding + reading frame mismatches). |
| `identity_score`      | Percentage identity computed as matches / (matches + mismatches).                |
| `annot1_start`        | Genomic start coordinate of the reference gene.                                  |
| `annot1_end`          | Genomic end coordinate of the reference gene.                                    |
| `annot2_start`        | Genomic start coordinate of the alternative gene.                                |
| `annot2_end`          | Genomic end coordinate of the alternative gene.                                  |
| `annot1_mRNA`         | Identifier of the selected reference mRNA used for the comparison.               |
| `annot2_mRNA`         | Identifier of the selected alternative mRNA used for the comparison.             |
| `CNC_mismatches_zones` | Genomic intervals where nucleotides are annotated as CDS in only one of the two annotations.       |
| `RF_mismatches_zones` | Genomic intervals where reading frames differ between annotations.               |
| `CNC_mismatches`       | Total length (in nucleotides) of coding/non-coding mismatches.                         |
| `RF_mismatches`       | Total length (in nucleotides) of reading frame mismatches.                       |
| `annot1_mRNA_number`  | Number of mRNAs annotated for the reference gene.                                |
| `annot2_mRNA_number`  | Number of mRNAs annotated for the alternative gene.                              |

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

| Column            | Description                                                                       |
| ----------------- | --------------------------------------------------------------------------------- |
| `Reference_locus` | Gene identifier in the reference annotation.                                      |
| `<alt>.locus`     | Identifier of the best matching gene in the corresponding alternative annotation. |
| `<alt>.identity`  | Identity score (%) for the corresponding gene pairing.                            |

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

[...]

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
