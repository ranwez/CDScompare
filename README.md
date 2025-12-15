
# CDScompR

This is the repository for the 2024 internship of Vetea Jacot in the GE²pop team of the research unit AGAP. The goal of the program CDScompR is to enable the comparison of two genome annotations in GFF format (one reference and one alternative) to determine the structure identity of the alternative to the reference. This program takes into account multiple caracteristics of the problem, for exemple the reading frames of a locus and overlapping loci. It creates two files: "log.txt" records all minor structure errors of the given input files, and "results.csv" lists the detailed results for each locus comparison. "CDScompR" stands for "CDS Comparison with Reading-frames".

This project also includes a script to enable the comparison of multiple alternative annotations to a single reference: CDSmulticompR. It relies on CDScompR for the comparison of each annotations pair, and is called in a similar way.



## Installation

### Requirements
One of the following must be available:
- Python ≥ 3.9 with `pip` 
- or Apptainer / Singularity

### Install from the Git repository with `pip` (users)

```
git clone git@github.com:ranwez/CDScompare.git
cd CDScompare
pip install .
```
*NB: Once installed, the repository is no longer required.*  

Successful installation can be tested with:
```
cdscompare --help
```

### Development installation (developers)

This method installs the package in editable mode and includes development dependencies (pytest).
```
git clone git@github.com:ranwez/CDScompare.git
cd CDScompare
pip install -e .[dev]
```

To run the tests:
```
python -m pytest -v
```

## Running the programs

To run CDScompR (one-to-one comparison), use the following command :

```
python path/to/CDScompR/script/CDScompR.py [ -h/--help -v/--verbose -o/--old_version -e/--exon_mode ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ]
```

where <reference_file_path> is the path to the file containing the reference annotation to analyze, and <alternative_file_path> is the path to the file containing the alternative annotation. These two parameters are required, while all others are not.

optional parameters :
 
-h/--help : displays script usage

-v/--verbose : displays messages indicating progress during execution

-o/--old_version : triggers use of the old program version which uses comparison of entire locus structure strings to compute identity, which can take up more memory during execution

-e/--exon_mode : triggers use of the exon coordinates instead of the CDS coordinates to compare the two annotations' loci

To run CDSmulticompR (multiple comparisons), use the same command with multiple calls to the '-a' option:

```
python path/to/CDScompR/script/CDSmulticompR.py [ -h/--help -v/--verbose -o/--old_version -e/--exon_mode ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> -a/--alternative <second_alternative_file_path> -a/--alternative <third_alternative_file_path> -a/--alternative ...] 
```


## Caution

The program expects files in GFF format (.gff or .gff3 extension) for which the structure lines are correctly ordered according to their relationships. It the structure lines order is not correct, an error is returned. If you need to compare disordered files, use the program "GFFcleaner" first to create new clean equivalent files: https://github.com/ranwez/GeneModelTransfer/blob/master/SCRIPT/VR/gff_cleaner.py (A new comprehensive annotation of leucine-rich repeat-containing receptors in rice, Gottin et al., 2021)


## Program test

A small example of reference and alternative public annotations are included in the repository and can be used to test and visualize program results before using it on real data. This example can be run using the command :

```
python path/to/CDScompR/script/CDScompR.py -r data/tests/program_test_ref.gff3 -a data/tests/program_test_alt.gff3
```


## Program's unit tests

An automatic unit testing script is included (data/tests/test_units.py) and can be used to check if an update didn't change the expected results of the program. Simply run :

```
pytest -vv
```

in the root (highest directory) of the repository. If the result of any function of the program is not what is expected, error messages should appear; if not, the program is working as intended. 

if pytest is not installed, run the following command :

```
pip install -U pytest
```
