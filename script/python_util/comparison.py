#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 09:38:04 2024

@author: vetea, ranwez
"""

import intervals_utils as iu
import pre_comparison as pc
from attrs import define, field
from typing import Optional
from locus import Locus

@define(slots=True, eq=True)
class MatchScore:
    """Class representing a score for a comparison."""
    genomic_overlap: int = field(default=0)
    identity: float = field(default=0.0)
    def is_better_than(self, other) -> bool:
        """Compare two scores based on genomic overlap and identity."""
        if self.genomic_overlap == 0:
            return False
        if other.genomic_overlap == 0:
            return True
        if self.identity == other.identity :
                return self.genomic_overlap > other.genomic_overlap
        else:
            return self.identity > other.identity

@define(slots=True, frozen=True, eq=True)
class MismatchInfo:
    """Class representing a mismatch zone between two annotations' loci."""
    zones: list = field(factory=list)  # List of mismatch zone coordinates with included bounds
    nb: int = field(default=0, init=False)  # Number of mismatches in the zone, calculated automatically
    
    def __attrs_post_init__(self):
        """Calculate number of mismatches after initialization."""
        if len(self.zones) % 2 != 0:
            raise ValueError("Mismatch zones must be defined by pairs of coordinates.")
            
        total = 0
        for i in range(0, len(self.zones), 2):
            total += self.zones[i+1] - self.zones[i] + 1
        
        object.__setattr__(self, 'nb', total)
    
    def __str__(self):
        """String representation of the mismatch information."""
        return f"Number of mismatches: {self.nb}, Zones: {self.zones}"
        
    @classmethod
    def create(cls, zone_bounds=None):
        """Factory method to create a MismatchInfo instance."""
        if zone_bounds is None:
            zone_bounds = []
        return cls(zones=zone_bounds.copy())

@define(slots=True, frozen=True, eq=True)
class MrnaMatchInfo:
    """Class representing the result of a comparison between two annotations' loci."""
     # Comparison metrics
    matches: int = field(default=0)
    mismatches_EI: MismatchInfo = field(factory=MismatchInfo)  # Exon-Intron mismatches
    mismatches_RF: MismatchInfo = field(factory=MismatchInfo)  # Exon-Intron mismatches
    genomic_overlap: int = field(default=0)  # Overlap in genomic coordinates
    ref_id: str = field(default='_') 
    alt_id: str =  field(default='_')
    score : MatchScore = field(factory=MatchScore, init=False)  # Score for the match
    
        
    def __attrs_post_init__(self):
        """Calculate identity after initialization."""
        total = self.matches + self.mismatches_EI.nb + self.mismatches_RF.nb
        if total == 0:
            object.__setattr__(self, 'score', MatchScore(self.genomic_overlap, 0))
        else:
            object.__setattr__(self, 'score', MatchScore(self.genomic_overlap, self.matches / total))

    def get_identity(self) -> float:
        return self.score.identity

    def has_better_identity_than(self, other: 'MrnaMatchInfo') -> bool:
        return self.score.is_better_than(other.score)



## reverses the coordinates of the mismatch zones returned by the functions
# compare_loci
#
# @see compare_loci()
#
# @param mismatch_zones Comparison mismatch zones as returned by compare_loci 
# (tuple of lists)
#
# @param cluster_end End coordinate of the last locus of the cluster
#
# @returns Returns a tuple of lists corresponding to the reversed zones
def reverse_coord(mismatch_zones, cluster_end):
    new_list_EI = []
    new_list_RF = []
    
    if mismatch_zones[0] != []:
        for bound in reversed(mismatch_zones[0]):
            new_list_EI.append(cluster_end-bound)
        
    if mismatch_zones[1] != []:
        for bound in reversed(mismatch_zones[1]):
            new_list_RF.append(cluster_end-bound)
            
    return (new_list_EI, new_list_RF)

def overlap(interval1, interval2) -> int:
        if interval1[1] < interval2[0] or interval2[1] < interval1[0]:
            return 0
        return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0]) + 1

## This function compares two annotations' loci returned by the function 
# get_gff_borders and creates for each pair of reference-alternative mRNAs 
# a comparison list detailing the identities and differences between the two
# annotations's codon position structure. It returns a MrnaMatchInfo instonce containing the 
# mRNA id (ref an alt) giving the highest identity and detail of their comparison
#
# @see get_gff_borders()
#
# @param ref_locus The reference annotation's locus (Locus class instance)
#
# @param alt_locus The alternative annotation's locus (Locus class instance)
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @returns Returns a MrnaMatchInfo 
# @see MrnaMatchInfo
#
# @remark This function doesn't expect any annotation to be a 'reference'
def compare_loci(ref_locus, alt_locus, verbose=False)-> MrnaMatchInfo:
    best_mRNA_comparison = MrnaMatchInfo()
    
    loci_overlap = overlap((ref_locus.start, ref_locus.end), (alt_locus.start, alt_locus.end)) 
    # if the loci are on different strands or don't overlap return 'null' values
    if ref_locus.direction != alt_locus.direction or loci_overlap == 0:
        return  best_mRNA_comparison
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
            matchInfo: MrnaMatchInfo = compute_matches_mismatches_EI_RF(mRNA_ref_id, mRNA_alt_id, mRNA_ref,  mRNA_alt, verbose)
            if matchInfo.has_better_identity_than(best_mRNA_comparison):
                best_mRNA_comparison = matchInfo
    
    return (best_mRNA_comparison)


## Computes the number of locus positions matching (both in CDS and same codon 
# position) or mismatching (one not in CDS or different codon position) in the
# given intervals
#
# @param mRNA_ref CDS position intervals of the reference (as a list)
#
# @param intervals_ref Instance of class 'OrderedIntervals' describing the 
# position intervals of the reference
#
# @param mRNA_alt CDS position intervals of the reference (as a list)
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see OrderedIntervals
#
# @returns Returns a MrnaMatchInfo providing the number of positions matching, mismatching
# for 'EI', mismatching for 'RF', the 'EI' mismatch zones, and the 'RF' 
# mismatch zones
#
# @remark 'EI' stands for 'Exon-Intron' and indicates one of the annotations
# is not in a CDS for the current position; 'RF' stands for 'Reading Frame' and
# indicates the two annotations are not in the same codon position at the 
# current position
#
# @remark This function was written by Vincent Ranwez
def compute_matches_mismatches_EI_RF(mRNA_ref_id:str, mRNA_alt_id:str, mRNA_ref, mRNA_alt, verbose)->MrnaMatchInfo:
    intervals_alt = iu.OrderedIntervals(mRNA_alt, True);
    intervals_ref = iu.OrderedIntervals(mRNA_ref, True);
    
    # get intersection and union of both intervals lists
    inter_mrna = intervals_ref.intersection(intervals_alt);
    union_mrna = intervals_ref.union(intervals_alt);
    
    # get exon and intron (EI) mismatches and mismatches zones (simple symmetric difference)
    mismatches_EI=union_mrna.total_length()-inter_mrna.total_length();
    diff_EI=intervals_ref.symmetric_difference(intervals_alt);
    
    # get intervals of intersection with their upper bounds
    inter_mrna_bounds=inter_mrna.get_intervals_with_included_ub();
    
    # get reading frames for each interval of the intersection
    inter_mrna_ref_rf=pc.get_reading_frame(mRNA_ref, inter_mrna_bounds, True)
    inter_mrna_alt_rf=pc.get_reading_frame(mRNA_alt, inter_mrna_bounds, True)
    
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
    mismatches_EI = MismatchInfo(diff_EI.get_intervals_with_included_ub())
    mismatches_RF = MismatchInfo(diff_RF)
    genomic_overlap = overlap((mRNA_ref[0], mRNA_ref[-1]), (mRNA_alt[0], mRNA_alt[-1]))
    return MrnaMatchInfo(matches, mismatches_EI, mismatches_RF, genomic_overlap, ref_id=mRNA_ref_id, alt_id=mRNA_alt_id)

def build_alignment_res(ref_locus, alt_locus, comparison, cluster, create_strings=False):
    """Build a result dictionary for an alignment match between two loci.
    
    Args:
        ref_locus: The reference locus object
        alt_locus: The alternative locus object
        comparison: MrnaMatchInfo object containing comparison results
        cluster_name: Name of the cluster (optional)
        create_strings: Whether to format mismatch zones as strings (for backward compatibility)
    
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
        num_mRNAs_ref = len(ref_locus.mRNAs.keys())
        result.update({
            "reference": ref_locus.name,
            "reference start": ref_locus.start,
            "reference end": ref_locus.end,
            "reference mRNA number": num_mRNAs_ref
        })
        if(comparison is None):
            def_ref_id =  list(ref_locus.mRNAs.keys())[0] if ref_locus.mRNAs else "_"
            result.update({
                "reference mRNA": def_ref_id 
            })
            
    if alt_locus is not None:
        num_mRNAs_alt = len(alt_locus.mRNAs.keys())
        result.update({
            "alternative": alt_locus.name,
            "alternative start": alt_locus.start,
            "alternative end": alt_locus.end,
            "alternative mRNA number": num_mRNAs_alt
        })
        if comparison is None:
            def_alt_id = list(alt_locus.mRNAs.keys())[0] if alt_locus.mRNAs else "_"
            result.update({
                "alternative mRNA": def_alt_id
            })

    if comparison is not None:
        # Extract values from the comparison
        identity = comparison.get_identity() 
        mismatch_zones = (comparison.mismatches_EI.zones, comparison.mismatches_RF.zones)
        
        # Handle reverse strand if needed
        if ref_locus.direction == "reverse":
            mismatch_zones = reverse_coord(mismatch_zones, cluster.get_end())

        # Convert mismatch zones to string representation for compatibility with old code
        if create_strings:
            mismatch_zones = "?" if mismatch_zones == ([], []) else str(mismatch_zones)
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

## This function compares the loci of the reference and alternative clusters
# returned by the construct_clusters function by assigning each locus to the 
# overlapping locus of the other annotation cluster which gives the highest
# computed identity. Each comparison results are written to a 'results' 
# dictionary detailing locus information.
#
# @see construct_clusters()
#
# @param cluster Instance of class 'Cluster' as returned by construct_clusters
#
# @param create_strings Boolean indicating wether to use the new compare_loci
# fucntion ('False') or the old_compare_loci function ('True')
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see compare_loci()
#
# @returns Returns a list of dictionaries detailing the locus information for 
# each locus comparison 
def annotation_match(cluster, create_strings=False, verbose=False):
    cluster_ref : list[Locus]  = cluster.get_loci()["ref"]
    cluster_alt : list [Locus] = cluster.get_loci()["alt"]
    cluster_name = cluster.name
    if verbose:
        print(f"matching annotations loci with each other for {cluster_name}")
    dyn_prog_matrix = [] # dynamic programmation matrix
    
    # if the loci stored in the cluster are on the reverse strand, reverse 
    # their mRNA's cds lists
    try:
        ref_is_rev = cluster_ref[0].direction == "reverse" 
    except IndexError:
        ref_is_rev = False
    try:
        alt_is_rev = cluster_alt[0].direction == "reverse" 
    except IndexError:
        alt_is_rev = False
    if ref_is_rev == True and alt_is_rev == True:
        cluster_end = cluster.get_end()
        for loc in cluster_ref:
            loc.reverse(cluster_end)
        for loc in cluster_alt:
            loc.reverse(cluster_end)
                 
    # Initialisation de la matrice de programmation dynamique avec des zéros
    dyn_prog_matrix = [[0 for j in range(len(cluster_alt) + 1)] for i in range(len(cluster_ref) + 1)]       
            
    # compute all internal values of the matrix by taking the maximum value of
    # the top 'cell', the left 'cell', and the top-left 'cell' summed with the
    # computed identity for its corresponding loci
    for i in range(1, len(cluster_ref)+1):
        for j in range(1, len(cluster_alt)+1):
            try:
                comparison = compare_loci(cluster_ref[i-1], cluster_alt[j-1], verbose)
            except Exception as e:
                ref_name = getattr(cluster_ref[i-1], 'name', '?')
                alt_name = getattr(cluster_alt[j-1], 'name', '?')
                print(f"Exception during loci comparison: ref_locus={ref_name}, alt_locus={alt_name}")
                print(f"Exception: {e}")
                raise
            
            dyn_prog_matrix[i][j] = max(dyn_prog_matrix[i-1][j],
                                        dyn_prog_matrix[i][j-1],
                                        dyn_prog_matrix[i-1][j-1] + comparison.get_identity())
        
    # retrieve best match alignment through backtracking
    i = len(cluster_ref)
    j = len(cluster_alt)
    results = []
    while i>0 and j>0:
        comparison = compare_loci(cluster_ref[i-1], cluster_alt[j-1], False)
        
        if (dyn_prog_matrix[i][j] == dyn_prog_matrix[i-1][j-1] + comparison.get_identity()): 
            results.append(build_alignment_res(cluster_ref[i-1], cluster_alt[j-1],comparison, cluster))
            i -= 1
            j -= 1
        elif( dyn_prog_matrix[i-1][j] == dyn_prog_matrix[i][j]):
            results.append(build_alignment_res(cluster_ref[i-1], None ,None, cluster))
            i -= 1
        else:
            results.append(build_alignment_res(None, cluster_alt[j-1],None, cluster))
            j -= 1

    # add loci remaining in the reference or alternative annotations
    while j > 0:
            results.append(build_alignment_res(None, cluster_alt[j-1],None, cluster))
            j -= 1
    while i > 0 :
            results.append(build_alignment_res(cluster_ref[i-1], None ,None, cluster))
            i -= 1
    # sort the results dictionary list so they are in order of locus start
    final_results = sorted(results, key=lambda d: d['reference'])
        
    return final_results
    
    
