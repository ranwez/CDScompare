#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:46 2024

@author: vetea
"""

import sys
import os
import re
import locus as lc


## This function retrieves and returns the id of the structure described from
# a line read from a GFF file.
#
# @param parsed_line The list of the line read from the file parsed through its
# tabulations (string 'split' method)
#
# @param delim Instance of class 'Pattern' ('re' library) compiled with the
# regular expression selecting all characters from the start of the line to
# the first delimiter character
#
# @param verbose If True, triggers display of more information messages.
# Default is 'False'
#
# @remark This function expects the file to be in GFF format
#
# @returns Returns the id of the structure described by the line
def get_structure_id(parsed_line, delim, verbose=False):

    # try to retrieve the last column, else return an error
    try:
        last_col = parsed_line[8]
    except IndexError:
        print("\nFile is not in GFF format (no ninth column)")
        sys.exit(1)

    # retrieve the id with the rest of the column text (commentaries)
    try:
        id_and_rest = last_col.split("=")[1]
    except IndexError:
        print("\nNo ID field found in last column (ninth column)")
        sys.exit(1)

    structure_id = delim.findall(id_and_rest)[0]

    return structure_id


## This function retrieves and returns the parent id of the structure
# described from a line read from a GFF file.
#
# @param parsed_line The list of the line read from the file parsed through its
# tabulations (string 'split' method)
#
# @param delim Instance of class 'Pattern' ('re' library) compiled with the
# regular expression selecting all characters from the start of the line to
# the first delimiter character
#
# @param verbose If True, triggers display of more information messages.
# Default is 'False'
#
# @remark This function expects the file to be in GFF format, and the
# parent id field to come after the id field
#
# @returns Returns the id of the parent of the structure described by the line
def get_parent_id(parsed_line, delim, verbose=False):

    # try to retrieve the last column, else return an error
    try:
        last_col = parsed_line[8]
    except IndexError:
        print("\nFile is not in GFF format (no ninth column)")
        sys.exit(1)

    # retrieve the id with the rest of the column text (commentaries)
    try:
        id_and_rest = last_col.split("=")[2]
    except IndexError:
        print("\nNo parent ID field found in last column (ninth column)")
        sys.exit(1)

    structure_id = ""

    structure_id = delim.findall(id_and_rest)[0]

    return structure_id


## This function expects a string corresponding to the file path of the GFF
# file to read, and returns a dictionary of dictionaries of instances of the
# class 'Locus', detailing all the relevant information for the gene and
# its mRNAs for each chromosome and DNA strand (direct/reverse)
#
# @see Locus
#
# @param path Path of the file to read
#
# @param verbose If True, triggers display of more information messages.
# Default is 'False'
#
# @param exon_mode Boolean indicating if the main comparison structures read
# from the file should be coding sequences (CDS, False) or exons (True).
# Default is 'False' (CDS comaprison)
#
# @remark If the parent ID of a CDS does not match the ID of the previous
# mRNA (indicating an incorrect file structure), an entry is added to a 'log'
# file but the function is not interrupted
#
# @return Returns a dictionary of dictionaries of instances of the class
# 'Locus', containing for each chromosome and DNA strand the information of
# the CDS borders of each mRNA of the gene, the start and end coordinates, the
# DNA strand on which the gene is predicted, and the locus ID
def get_gff_borders(path, verbose=False, exon_mode=False):

    # specifiy the desired main structure to be read
    if exon_mode:
        exon_or_cds = "exon"
    else:
        exon_or_cds = "CDS"

    try:
        file = open(path, "r") # the file to read
    except FileNotFoundError:
        print(f"\nget_gff_borders() function error : file '{path}' does not exist or cannot be read from\n")
        sys.exit(2)

    # try to open the structure errors log file
    try:
        log = open("./results/log.txt", "w")
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        log = open("./results/log.txt", "w")

    all_loci = {} # return dictionary
    locus_id = "" # the current gene being analyzed
    mRNA_id = "" # the current mRNA being analyzed
    line_index = 1 # number of the file line currently read
    locus = lc.Locus() # Initialisation of the Locus class instance to construct
    loci = None

    delim = re.compile("^[^,?;:/!$%@#~&\n\t]+")

    for line in file:

        parsed_line = line.split("\t")

        if parsed_line[0] == "contig":
            continue

        if str(parsed_line[2]) == "gene":

            # if there was a previous gene, add it to return dictionary and
            # delete it. Then we create a new Locus instance
            # (if the gene is on the reverse strand, reverse its border list
            # before adding it to the dictionary)
            if loci != None:

                # if there was a previous gene, but its borders list is empty,
                # return an error
                if locus_id != "" and locus.mRNAs[mRNA_id] == []:
                    print(f"\nLine {line_index} = get_gff_borders() function error : no coding sequence (CDS) could be found for the previous mRNA '{mRNA_id}' in file {path}\nit is possible the file has an incorrect features order. You can clean it using https://github.com/ranwez/GeneModelTransfer/blob/master/SCRIPT/VR/gff_cleaner.py\n")
                    sys.exit(1)

                loci[locus_id] = locus
            locus = lc.Locus() # new instance for the locus being read
            mRNA_id = "" # identifier of the mRNA currently being read
            parent_id = "" # identifier of the parent of the current structure
            locus_id = "" # identifier of the current locus

        chrm = parsed_line[0] # name of the chromosome
        if chrm == "contig":
            continue
        strand = "direct" if parsed_line[6] == "+" else "reverse"
        try:
            loci = all_loci[chrm+"_"+strand]
        except KeyError:
            all_loci[chrm+"_"+strand] = {}
            loci = all_loci[chrm+"_"+strand]

        # if we encounter a new gene which is not from a contig,
        # we get its information and create an instance in 'loci'
        if str(parsed_line[2]) == "gene":

            locus.direction = strand

            # we retrieve the id of the locus
            locus_id = get_structure_id(parsed_line, delim, verbose)
            locus.name = locus_id

            # we retrieve the start and end coordinates of the locus
            start = int(parsed_line[3])
            locus.start = start # start coordinate of the locus
            end = int(parsed_line[4])
            locus.end = end # end coordinate of the locus

            if verbose :
                print("Reading the locus " + locus_id)

        # if we encounter a new mRNA, we get its ID and create a key
        # in the instance's mRNAs attribute with an empty list
        elif str(parsed_line[2]) == "mRNA":

            mRNA_id = get_structure_id(parsed_line, delim, verbose)
            locus.mRNAs[mRNA_id] = [] # list of CDS coordinates of the mRNA

        # if we encounter a CDS, we add its start and end positions to the
        # corresponding mRNA key in the locus' mRNAs attribute
        elif str(parsed_line[2]) == exon_or_cds:
            parent_id = get_parent_id(parsed_line, delim, verbose)

            if parent_id != mRNA_id:
                print(f"\nIncorrect file structure (Parent of CDS is not previous mRNA) in file {path}. See 'log.txt' for more information")
                log.write("Line " + str(line_index) + " : CDS parent ID (" + parent_id + ") does not match last mRNA ID (" + locus_id +")\n")

            if mRNA_id == '':
                print(f"\nLine {line_index} = get_gff_borders() function error : CDS has been found before any mRNA in file {path}\nit is possible the file has an incorrect features order. You can clean it using https://github.com/ranwez/GeneModelTransfer/blob/master/SCRIPT/VR/gff_cleaner.py\n")
                sys.exit(1)

            locus.mRNAs[mRNA_id].append(int(parsed_line[3]))
            locus.mRNAs[mRNA_id].append(int(parsed_line[4]))

        line_index += 1 # increment the line indicator (for debug purposes)

    if locus_id == "": # if the locus_id still has default value, return error
        print(f"\nget_gff_borders() function error : no locus could be found in file '{path}'\n")
        sys.exit(1)

    file.close()
    log.close()
    loci[locus_id] = locus # add the last locus to the return dictionary

    # we return the entire dictionary with all loci (of all chromosomes)
    return all_loci
