#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def format_mismatch_zones(zones):
    """Formate une liste de coordonnées de mismatch en chaîne lisible"""
    return " ".join(f"[{zones[i]}//{zones[i+1]}]" for i in range(0, len(zones), 2))

def write_results(all_results, alt_name, out_dir):
    # gene stats [gene found in both; found only in ref; found only in alt]
    full_annotation_stat = [0, 0, 0]

    # try to open the results file (named after the alternative file name)
    filename=f"{out_dir}/{alt_name}.csv"
    try:
        results_file = open(filename, "w")
    except FileNotFoundError:
        os.mkdir(out_dir) # create 'results' subdirectory
        results_file = open(filename, "w")

    results_file.write("Chromosome,Cluster name,Reference locus,Alternative locus,Comparison matches,Comparison mismatches,Identity score (%),Reference start,Reference end,Alternative start,Alternative end,Reference mRNA,Alternative mRNA,Exon_intron (EI) non-correspondance zones,Reading frame (RF) non-correspondance zones,Exon_Intron (EI) mismatches,Reading Frame (RF) mismatches,reference mRNA number,alternative mRNA number\n")

    for dna_mol, results in all_results.items():
        chromosome_strand_stat = [0,0,0]

        for cluster in results:
            for loc in cluster:
                if loc['mismatch zones'] not in ["_", "?"]:
                    # convert mismatch zones so commas don't modify the CSV output
                    mismatch_EI = format_mismatch_zones(loc['mismatch zones'][0])
                    mismatch_RF = format_mismatch_zones(loc['mismatch zones'][1])

                else:
                    mismatch_EI = loc['mismatch zones']
                    mismatch_RF = loc['mismatch zones']

                if loc['mismatch/match'] == []:
                    results_file.write(f"{dna_mol},{loc['cluster name']},{loc['reference']},{loc['alternative']},_,_,{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']},{loc['reference mRNA']},{loc['alternative mRNA']},_,_,_,_,{loc['reference mRNA number']},{loc['alternative mRNA number']}\n")
                    if loc['reference'] == '~':
                        chromosome_strand_stat[2] += 1
                    else:
                        chromosome_strand_stat[1] += 1
                else:
                    results_file.write(f"{dna_mol},{loc['cluster name']},{loc['reference']},{loc['alternative']},{loc['mismatch/match'][0]},{loc['mismatch/match'][1]+loc['mismatch/match'][2]},{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']},{loc['reference mRNA']},{loc['alternative mRNA']},{mismatch_EI},{mismatch_RF},{loc['mismatch/match'][1]},{loc['mismatch/match'][2]},{loc['reference mRNA number']},{loc['alternative mRNA number']}\n")
                    chromosome_strand_stat[0] += 1

        print(f"\nNumber of loci of {dna_mol}:\n- found in both annotations : {chromosome_strand_stat[0]}\n- found only in the reference : {chromosome_strand_stat[1]}\n- found only in the alternative : {chromosome_strand_stat[2]}\n")
        full_annotation_stat = [sum(x) for x in zip(full_annotation_stat, chromosome_strand_stat)]

    results_file.close()

    results_file = open(f"{out_dir}/{alt_name}.txt", "w")
    results_file.write(f"\nNumber of loci (whole data):\n- found in both annotations : {full_annotation_stat[0]}\n- found only in the reference : {full_annotation_stat[1]}\n- found only in the alternative : {full_annotation_stat[2]}\n")
    results_file.close()
    print(f"\nNumber of loci (whole data):\n- found in both annotations : {full_annotation_stat[0]}\n- found only in the reference : {full_annotation_stat[1]}\n- found only in the alternative : {full_annotation_stat[2]}\n")


## Writes the results returned by the function multicomp into a results 
# synthesis CSV file detailing the loci identity for each alternative
#
# @see multicomp()
#
# @param multi_results List of loci identities (dictionaries), as returned 
# by multicomp
#
# @param ref_path Path to the reference annotation GFF file
def write_multi_results(multi_results, ref_path, alt_paths):
    
    ref_name = os.path.basename(ref_path).split(".")[0]
    
    # get each reference key of the first dictionary (since the reference is 
    # the same for all comparisons) and use it to retrieve the corresponding
    # comparisons for all result dictionaries, then write them in a synthesis
    # CSV file
    
    # try to open the results file
    try:
        results_file = open(f"./results/synthesis_{ref_name}.csv", "w")
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        results_file = open(f"./results/synthesis_{ref_name}.csv", "w")
            
    # write the header    
    
    header = "Reference_locus"
    
    for alt in alt_paths:
        alt_name =  os.path.basename(alt).split(".")[0]
        header += f",{alt_name} locus,{alt_name} identity"
        
    results_file.write(header+"\n")
    
    # write the results
    ref_keys = set();
    for result in multi_results:
        ref_keys=ref_keys.union(result.keys())
    
    for ref_key in sorted(ref_keys):
        line = ref_key
        for alt in multi_results:
            if ref_key in alt:
                line += f",{alt[ref_key][0]},{alt[ref_key][1]}"
            else:
                line += ",~,0.0"
        results_file.write(line+"\n")
                    
    results_file.close()