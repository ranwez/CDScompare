# -*- coding: utf-8 -*-

from pathlib import Path
from script.python_util.annotation import AnnotationSet


def format_mismatch_zones(zones):
    """Format a list of mismatch coordinates into a readable string"""
    return " ".join(f"[{zones[i]}//{zones[i+1]}]" for i in range(0, len(zones), 2))

def write_results(all_results: dict, csv_path: Path, txt_path: Path):
    # gene stats [gene found in both; found only in ref; found only in alt]
    full_annotation_stat = [0, 0, 0]

    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with open(csv_path, "w") as results_file:
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

    with open(txt_path, "w") as results_file:
        results_file.write(f"\nNumber of loci (whole data):\n- found in both annotations : {full_annotation_stat[0]}\n- found only in the reference : {full_annotation_stat[1]}\n- found only in the alternative : {full_annotation_stat[2]}\n")

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
def write_multi_results(multi_results: list[dict], annotations: AnnotationSet, out_dir: Path):
    
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = annotations.synthesis_filename(out_dir)

    # get each reference key of the first dictionary (since the reference is 
    # the same for all comparisons) and use it to retrieve the corresponding
    # comparisons for all result dictionaries, then write them in a synthesis
    # CSV file
    
    with open(csv_path, "w") as results_file:        
        # write the header    
        header = "Reference_locus"
        for alt in annotations.alts:
            header += f",{alt.id} locus,{alt.id} identity"
        results_file.write(header+"\n")
        
        # write the results
        ref_keys = set()
        for result in multi_results:
            ref_keys=ref_keys.union(result.keys())
        
        for ref_key in sorted(ref_keys):
            line = ref_key
            for alt_result in multi_results:
                if ref_key in alt_result:
                    line += f",{alt_result[ref_key][0]},{alt_result[ref_key][1]}"
                else:
                    line += ",~,0.0"
            results_file.write(line+"\n")
                    
