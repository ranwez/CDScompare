#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 09:38:04 2024

@author: vetea, ranwez
"""

from locus import Locus
from cluster import Cluster
from match import MatchScore, MrnaMatchInfo, MismatchInfo

def get_reading_frame(cds_bounds, area_bounds, phase_first_CDS=0, verbose=False):
    nb_nt = (3 - phase_first_CDS) % 3
    nb_nt_in_cds=0;
    cdsb = 0
    cds_len = len(cds_bounds)
    reading_frames =[]
    for area_start in area_bounds[::2]:
        # Move to CDS that includes current area start
        while cdsb + 1 < cds_len and area_start > cds_bounds[cdsb + 1]:
            nb_nt += cds_bounds[cdsb + 1] - cds_bounds[cdsb] + 1
            cdsb += 2
        # Now set the reading frame of the current area
        nb_nt_in_cds = area_start - cds_bounds[cdsb] + 1
        reading_frames.append((nb_nt + nb_nt_in_cds - 1) % 3 + 1)
    return reading_frames


def reverse_coord(mismatch_zones, cluster_end):
    """Reverse the coordinates of mismatch zones relative to the cluster end."""
    return (
        [cluster_end - b for b in reversed(mismatch_zones[0])] if mismatch_zones[0] else [],
        [cluster_end - b for b in reversed(mismatch_zones[1])] if mismatch_zones[1] else []
    )

def overlap(interval1: tuple[int, int], interval2: tuple[int, int]) -> int:
    """Calculate the overlap between two intervals [start, end]."""
    if interval1[1] < interval2[0] or interval2[1] < interval1[0]:
        return 0
    return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0]) + 1


def compare_loci(ref_locus: Locus, alt_locus: Locus) -> MrnaMatchInfo:
    """Compare two loci and return the best mRNA comparison."""
    best_mRNA_comparison = MrnaMatchInfo()
    if overlap((ref_locus.start, ref_locus.end), (alt_locus.start, alt_locus.end)) == 0:
        return best_mRNA_comparison

    nb_alt_mrnas = len(alt_locus.mRNAs)
    for i_ref in range(len(ref_locus.mRNAs)):
        for i_alt in range (nb_alt_mrnas):
            matchInfo: MrnaMatchInfo = compute_matches_mismatches_EI_RF(i_ref, i_alt, ref_locus, alt_locus)
            if matchInfo.has_better_identity_than(best_mRNA_comparison):
                best_mRNA_comparison = matchInfo
    return best_mRNA_comparison

def compute_matches_mismatches_EI_RF(i_ref, i_alt, ref:Locus, alt:Locus)->MrnaMatchInfo:
    """Compute MrnaMatchInfo between two mRNAs."""
    intervals_ref= ref.mrna_intervals[i_ref]
    intervals_alt = alt.mrna_intervals[i_alt]
    mrna_ref= ref.mRNAs[i_ref]
    mrna_alt= alt.mRNAs[i_alt]

    # get intersection and union and symDiff of both intervals lists
    inter_mrna, union_mrna, diff_EI = intervals_ref.inter_union_symdiff(intervals_alt)

    # get exon and intron (EI) mismatches and mismatches zones (simple symmetric difference)
    mismatches_EI=union_mrna.total_length()-inter_mrna.total_length()

    inter_mrna_bounds=inter_mrna.as_list_with_included_ub()

    inter_mrna_ref_rf=get_reading_frame(mrna_ref, inter_mrna_bounds, ref.phases[i_ref], True)
    inter_mrna_alt_rf=get_reading_frame(mrna_alt, inter_mrna_bounds, alt.phases[i_alt], True)

    # for each intersection interval of the reference and alternative, compare
    # the reading frames to determine match or RF mismatch
    matches : int=0 
    mismatches_RF : int = 0
    diff_RF : list = [] # reading frame (RF) mismatches zones
    for interval_id in range(0, len(inter_mrna_ref_rf)):
        interval_lg = inter_mrna_bounds[2*interval_id+1] - inter_mrna_bounds[2*interval_id] + 1
        if inter_mrna_ref_rf[interval_id] != inter_mrna_alt_rf[interval_id]:
            mismatches_RF += interval_lg
            diff_RF.append(inter_mrna_bounds[2*interval_id])
            diff_RF.append(inter_mrna_bounds[2*interval_id+1])
        else:
            matches += interval_lg
    mismatches_EI = MismatchInfo(diff_EI.as_list_with_included_ub())
    mismatches_RF = MismatchInfo(diff_RF)
    genomic_overlap = overlap((mrna_ref[0], mrna_ref[-1]), (mrna_alt[0], mrna_alt[-1]))
    return MrnaMatchInfo(matches, mismatches_EI, mismatches_RF, genomic_overlap, ref_id=ref.ids[i_ref], alt_id=alt.ids[i_alt])

def build_alignment_res(ref_locus, alt_locus, comparison, cluster, reversed):
    """Build a result dictionary for an alignment match between two loci.
    
    Args:
        ref_locus: The reference locus object
        alt_locus: The alternative locus object
        comparison: MrnaMatchInfo object containing comparison results
        cluster_name: Name of the cluster (optional)
        direction: Direction of the locus ('direct' or 'reverse')
    
    Returns:
        Dictionary containing all alignment information
    """
    
    #build default res and add available information
    result= {"reference" : "~",
            "reference start" : "_",
            "reference end" : "_",
            "alternative" : "~",
            "alternative start" : "_",
            "alternative end" : "_",
            "reference mRNA" : "_",
            "alternative mRNA" : "_",
            "mismatch/match" : [],
            "identity" : 0.0,
            "mismatch zones" : "_",
            "cluster name" : cluster.name,
            "reference mRNA number" : "_",
            "alternative mRNA number" : "_"}
    
    if ref_locus is not None:
        num_mRNAs_ref = len(ref_locus.mRNAs)
        result.update({
            "reference": ref_locus.name,
            "reference start": ref_locus.start,
            "reference end": ref_locus.end,
            "reference mRNA number": num_mRNAs_ref
        })
        if(comparison is None):
            def_ref_id = ref_locus.ids[0] if ref_locus.ids else "_"
            result.update({
                "reference mRNA": def_ref_id 
            })
            
    if alt_locus is not None:
        num_mRNAs_alt = len(alt_locus.mRNAs)
        result.update({
            "alternative": alt_locus.name,
            "alternative start": alt_locus.start,
            "alternative end": alt_locus.end,
            "alternative mRNA number": num_mRNAs_alt
        })
        if comparison is None:
            def_alt_id = alt_locus.ids[0] if alt_locus.ids else "_"
            result.update({
                "alternative mRNA": def_alt_id
            })

    if comparison is not None:
        # Extract values from the comparison
        identity = comparison.get_identity() 
        mismatch_zones = (comparison.mismatches_EI.zones, comparison.mismatches_RF.zones)
        
        # Handle reverse strand if needed
        if reversed :
            mismatch_zones = reverse_coord(mismatch_zones, cluster.end)

        result.update({
            "reference mRNA": comparison.ref_id,
            "alternative mRNA": comparison.alt_id,
            "mismatch/match": [comparison.matches, 
                               comparison.mismatches_EI.nb, 
                               comparison.mismatches_RF.nb],
            "identity": round(identity * 100, 1) if identity is not None else 0.0,
            "mismatch zones": mismatch_zones
        })
    return result

def annotation_match(cluster: Cluster, reversed:bool, mode_align: bool):
    cluster_ref : list[Locus]  = cluster.loci_ref
    cluster_alt : list [Locus] = cluster.loci_alt
    results = []
    # handle trivial cases
    if not cluster_ref:
        for alt_locus in cluster_alt:
            results.append(build_alignment_res(None, alt_locus, None, cluster, reversed))
        return results
    if not cluster_alt:
        for ref_locus in cluster_ref:
            results.append(build_alignment_res(ref_locus, None, None, cluster, reversed))
        return results
    if reversed:
        cluster.reverse_loci_coord()  
    cluster.set_intervals()  
    if len(cluster_ref) == 1 and len(cluster_alt) == 1:
        comparison = compare_loci(cluster_ref[0], cluster_alt[0])
        if(comparison.score.genomic_overlap > 0):
            results.append(build_alignment_res(cluster_ref[0], cluster_alt[0],comparison, cluster, reversed))
        else: # if best is (0,0) dont display a match
            results.append(build_alignment_res(cluster_ref[0], None ,None, cluster,reversed))
            results.append(build_alignment_res(None, cluster_alt[0],None, cluster,reversed))
        return results
    
    if(mode_align):
        return annot_match_alignment(cluster, reversed)
    else:
        return annot_match_all(cluster, reversed)

def annot_match_all(cluster: Cluster, reversed:bool):
    cluster_ref : list[Locus]  = cluster.loci_ref
    cluster_alt : list [Locus] = cluster.loci_alt
    results = []
    paired_ref= set()
    paired_alt= set()
    # handle all pairs of loci in the cluster
    for i_ref in range(len(cluster_ref)):
        for i_alt in range(len(cluster_alt)):
            comparison = compare_loci(cluster_ref[i_ref], cluster_alt[i_alt])
            if(comparison.score.genomic_overlap > 0):
                results.append(build_alignment_res(cluster_ref[i_ref], cluster_alt[i_alt],comparison, cluster, reversed))
                paired_alt.add(i_alt)
                paired_ref.add(i_ref)
    for i_ref in range(len(cluster_ref)):
        if i_ref not in paired_ref:
            results.append(build_alignment_res(cluster_ref[i_ref], None ,None, cluster,reversed))
        for i_alt in range(len(cluster_alt)):
            if i_alt not in paired_alt:
                results.append(build_alignment_res(None, cluster_alt[i_alt],None, cluster,reversed))
    return results 


def annot_match_alignment(cluster: Cluster, reversed:bool):
    cluster_ref : list[Locus]  = cluster.loci_ref
    cluster_alt : list [Locus] = cluster.loci_alt
    results = []

    # handle general case using dynamic programming 
    dyn_prog_matrix = [[MatchScore() for j in range(len(cluster_alt) + 1)] for i in range(len(cluster_ref) + 1)]       
            
    for i in range(1, len(cluster_ref)+1):
        for j in range(1, len(cluster_alt)+1):
            comparison : MrnaMatchInfo = compare_loci(cluster_ref[i-1], cluster_alt[j-1])
            dyn_prog_matrix[i][j] = MatchScore.max3(dyn_prog_matrix[i-1][j],
                                        dyn_prog_matrix[i][j-1],
                                        MatchScore.add(dyn_prog_matrix[i-1][j-1], comparison.score))
        
    # retrieve best match alignment through backtracking
    i = len(cluster_ref)
    j = len(cluster_alt)
    
    while i>0 and j>0:
        comparison = compare_loci(cluster_ref[i-1], cluster_alt[j-1])
        
        if (dyn_prog_matrix[i][j] == MatchScore.add(dyn_prog_matrix[i-1][j-1], comparison.score)):
            if(comparison.score.genomic_overlap > 0):
                results.append(build_alignment_res(cluster_ref[i-1], cluster_alt[j-1],comparison, cluster, reversed))
            else: # if best is (0,0) dont display a match
                results.append(build_alignment_res(cluster_ref[i-1], None ,None, cluster,reversed))
                results.append(build_alignment_res(None, cluster_alt[j-1],None, cluster,reversed))
            i -= 1
            j -= 1
        elif( dyn_prog_matrix[i-1][j] == dyn_prog_matrix[i][j]):
            results.append(build_alignment_res(cluster_ref[i-1], None ,None, cluster,reversed))
            i -= 1
        else:
            results.append(build_alignment_res(None, cluster_alt[j-1],None, cluster,reversed))
            j -= 1

    while j > 0:
            results.append(build_alignment_res(None, cluster_alt[j-1],None, cluster,reversed))
            j -= 1
    while i > 0 :
            results.append(build_alignment_res(cluster_ref[i-1], None ,None, cluster,reversed))
            i -= 1
   
    # sort the results dictionary list so they are in order of locus start, inversing order is enough
    #final_results = sorted(results, key=lambda d: d['reference'])
    final_results = results[::-1]  
    return final_results
    
    
