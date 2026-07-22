# CDScompare

**CDScompare performs the quantitative comparison of structural genome annotations (GFF3) of the same genome.**  
For each pair of overlapping genes, it computes a frame-aware CDS similarity score based on coding-region agreement and reading-frame conservation.  
It is implemented in Python and distributed as a command-line tool.

The tool supports:

- pairwise comparison between two genome annotations
- multi-comparison of several annotations against a common reference
- two pairing strategies (best or all)

## Installation

### Requirements

One of the following must be available:

- Python ≥ 3.10 with `pip` or `pipx`
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

In particular, CDScompare assumes that all input annotation files:

* describe the **same genome assembly**;
* follow standard **GFF3 conventions**, including:
  * tab-delimited records with 9 columns;
  * `gene`, `mRNA` and `CDS` features organized in a consistent `gene -> mRNA -> CDS` hierarchy;
  * valid `ID` and `Parent` attributes linking genes, transcripts and CDS features;
  * unique gene and mRNA identifiers within each file;
  * coherent genomic coordinates;
  * non-overlapping CDS features within each transcript;
  * valid CDS phases (`0`, `1`, `2`, or `.`);
* use **compatible seqids** in column 1 across annotations.

### Seqid compatibility

The seqids in column 1 of the two GFF files must be compatible. CDScompare can harmonize an optional leading `chr` prefix when this improves the number of shared seqids and does not merge distinct seqids within the same annotation. For example, annotations using `chr1`, `chr2`, `chr3` can therefore be compared to annotations using `1`, `2`, `3`.

More complex naming differences, such as `NC_000001.11` vs `chr1`, `chromosome_1` vs `1`, or scaffold-specific naming schemes, must be harmonized by the user before running CDScompare.

### Recommended preprocessing with AGAT

Before running CDScompare, we strongly recommend standardizing each input annotation with [AGAT](https://github.com/NBISweden/AGAT) (Dainat, 2022). For example:

```bash
agat_convert_sp_gxf2gxf.pl -g annotation.gff3 -o annotation.agat.gff3
```

AGAT can help standardize GFF/GTF files, sort features, fix duplicated IDs, add missing `ID` and `Parent` attributes, add missing parent features when possible, and produce a more consistent GFF3 structure. This preprocessing step is especially recommended when input files come from different annotation tools or databases.

### Scope and assumptions

- Gene overlap is detected based on **gene genomic coordinates**.
- Similarity scores are computed **only from CDS features**.
- Only gene models represented as a `gene -> mRNA -> CDS` hierarchy are considered. Genes without CDS-containing mRNAs, transcripts without CDS features, and features not represented as `gene` entries are ignored.
- Alternative splicing is handled by selecting the **best-matching mRNA pair** for each gene pairing.

---

## Command-line interface

```bash
cdscompare [OPTIONS] ANNOT1_GFF ANNOT2_GFF [ANNOT3_GFF ...]
```

### Positional arguments

| Argument         | Description                                                                            |
| ---------------- | -------------------------------------------------------------------------------------- |
| `ANNOT1_GFF`     | First annotation file, reported as annotation 1 in pairwise output                     |
| `ANNOT2_GFF`     | Second annotation file, reported as annotation 2 in pairwise output                    |
| `ANNOT3_GFF ...` | Optional additional annotation files, each compared independently against `ANNOT1_GFF` |

> **Note:**
>
> - At least **two GFF files** are required.
> - In pairwise comparisons (exactly two GFF files), similarity scores and gene pairings are invariant to the input file order.
> - GFF input file basenames must be unique (used as annotation identifiers).

### Options

| Option               | Description                                                           |
| -------------------- | --------------------------------------------------------------------- |
| `-d, --out_dir`      | Output directory where result files are written (default: `results`). |
| `-p, --pairing_mode` | Pairing strategy used within clusters of overlapping genes. Possible values are:<br>• `best` (default): selects a globally optimal gene pairing using dynamic programming.<br>• `all`: reports all overlapping gene pairings without global optimization. |

## Output files

### Pairwise comparison

For each pairwise comparison, CDScompare produces two output files.

#### 1. Detailed comparison file (`<annotation1>_vs_<annotation2>.csv`)

This file contains one row per reported gene pair or unpaired gene, with the following columns:

| Column                 | Description                                                                                      |
| ---------------------- | ------------------------------------------------------------------------------------------------ |
| `seqid_strand`         | Input GFF seqid followed by the `_direct` or `_reverse` strand suffix                            |
| `cluster`              | Identifier of the cluster of overlapping genes                                                   |
| `annot1_gene`          | Gene identifier in annotation 1                                                                  |
| `annot2_gene`          | Gene identifier in annotation 2                                                                  |
| `matches`              | Number of nucleotide positions that are coding in both annotations and in the same reading frame |
| `mismatches`           | Total number of mismatched nucleotides (coding/non-coding + reading frame mismatches)            |
| `similarity_score`     | Frame-aware CDS similarity score (%), computed as 100 × matches / (matches + mismatches)         |
| `annot1_start`         | Start coordinate of the annotation 1 gene                                                        |
| `annot1_end`           | End coordinate of the annotation 1 gene                                                          |
| `annot2_start`         | Start coordinate of the annotation 2 gene                                                        |
| `annot2_end`           | End coordinate of the annotation 2 gene                                                          |
| `annot1_mRNA`          | Identifier of the selected mRNA from annotation 1 used for the comparison                        |
| `annot2_mRNA`          | Identifier of the selected mRNA from annotation 2 used for the comparison                        |
| `C_NC_mismatch_zones`  | Genomic intervals where nucleotides are annotated as CDS in only one of the two annotations      |
| `RF_mismatch_zones`    | Genomic intervals where reading frames differ between annotations                                |
| `C_NC_mismatches`      | Total length (in nucleotides) of coding/non-coding mismatches                                    |
| `RF_mismatches`        | Total length (in nucleotides) of reading frame mismatches                                        |
| `annot1_mRNA_count`    | Number of CDS-containing mRNAs for the corresponding gene in annotation 1                        |
| `annot2_mRNA_count`    | Number of CDS-containing mRNAs for the corresponding gene in annotation 2                        |

Special values:

- `_` : undefined or not applicable
- `~` : no corresponding gene reported for this comparison. This occurs when:
  - no overlapping gene exists in the other annotation, or
  - the gene was paired with a different gene in `best` pairing mode.

#### 2. Summary file (`<annotation1>_vs_<annotation2>.txt`)

This file contains overall counts for the corresponding pairwise comparison:

```
Comparison summary for <annotation1> vs <annotation2>:
- reported gene pairs: X
- unpaired genes in annotation 1: Y
- unpaired genes in annotation 2: Z
```

### Multi-comparison mode

When more than two annotation files are provided, the first annotation is used as a common reference and compared independently with each subsequent annotation.

For each comparison, CDScompare generates the same pairwise CSV and summary files described above. In these files:

- annotation 1 refers to the first input annotation;
- annotation 2 refers to the annotation currently compared with it.

An additional global summary file (`synthesis_<annotation1>.csv`) is generated, with one line per gene from the reference annotation, and two columns for each compared annotation:

| Column                           | Description                                                 |
| -------------------------------- | ----------------------------------------------------------- |
| `<reference>_gene`                  | Gene identifier in annotation 1                             |
| `<annotation>_gene`              | Best-matching gene from the compared annotation             |
| `<annotation>_similarity_score`  | Similarity score (%) for the corresponding gene pairing     |

## Citation

CDScompare has been accepted for publication in Plant Methods. Full citation details will be added when available.

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
