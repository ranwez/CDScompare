
# Annotation Comparison

This is the repository for the 2024 internship of Vetea Jacot in the GE²pop team of the research unit AGAP. The goal of the program defined in this repo is to enable the comparison of two genome annotations in GFF format (one reference and one alternative) to determine the structure identity of the alternative to the reference. This program takes into account multiple caracteristics of the problem, for exemple the reading frames of a locus and overlapping loci. It creates two files: "log.txt" records all minor structure errors of the given input files, and "results.csv" lists the detailed results for each locus comparison.


## Running the program

Run following command :

```
path/to/main.py [ -h/--help -d/--debug -v/--verbose -o/--old_version ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ]
```

where <reference_file_path> is the path to the file containing the reference annotation to analyze, and <alternative_file_path> is the path to the file containing the alternative annotation. These two parameters are required, while all others are not.

options:
 
-h/--help : displays script usage

-d/--debug : displays debug information during execution

-v/--verbose : displays messages indicating progress during execution

-o/--old_version : triggers use of the old program version which uses comparison of entire locus structure strings to compute identity, which can take up more memory during execution


## Testing the program

An automatic unit testing script is included (data/tests/test_units.py) and can be used to check if an update didn't change the expected results of the program. Simply run :

```
pytest
```

in the root (highest directory) of the repository. If the result of any function of the program is not what is expected, error messages should appear; if not, the program is working as intended. 
