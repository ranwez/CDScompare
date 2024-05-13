
# Annotation Comparison

This is the repository for the 2024 internship of Vetea Jacot in the GE²pop team of the research unit AGAP. The goal of the program defined in this repo is to enable the comparison of multiple genome annotations to determine the closest automatic annotation to a reference. This program takes into account multiple caracteristics of the problem, for exemple the reading frames of a locus.

## Running the program

Run following command :

```
path/to/main.py [ -h/--help -d/--debug -v/--verbose ] [ -i/--input <input_folder_path> ] [ -r/--reference <reference_file_name> ]
```

where <input_folder_path> is the path to the folder containing all GFF files to analyze (including the reference file), and <reference_file_name> is the name of the reference file *with its extension* (exemple : 'basic_test'). These two parameters are required, while all others are not.

## Testing the program

An automatic unit testing script is included (data/tests/test_basic.py) and can be used to check if an update didn't change the expected results of the program. Simply run :

```
pytest
```

in the root (highest directory) of the repository. If the result of any function of the program is not what is expected, error messages should appear; if not, the program is working as intended. 
