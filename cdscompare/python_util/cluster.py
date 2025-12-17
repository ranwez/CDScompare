# -*- coding: utf-8 -*-

from cdscompare.python_util.locus import Locus

### Class representing a cluster of overlapping loci
class Cluster:
    __slots__ = ('name', 'loci_ref', 'loci_alt', 'end')
    def __init__(self, name, loci_ref=[], loci_alt=[], end=-1):
        self.name = name
        self.loci_ref = loci_ref
        self.loci_alt = loci_alt
        self.end = end
        
    def reverse_loci_coord(self):
        """
        Reverse the coordinates of all loci in the cluster.
        This is used for loci on the reverse strand.
        """
        for locus in self.loci_ref:
            locus.reverse(self.end)
        for locus in self.loci_alt:
            locus.reverse(self.end)

    def set_intervals(self) :
        """Pre compute the mRNA intervals for all loci in the cluster."""
        for locus in self.loci_ref:
            locus.set_mrna_intervals()
        for locus in self.loci_alt:
            locus.set_mrna_intervals()

          
def build_cluster_list_from_Locus(ref_loci, alt_loci, dna_mol):
    """
    Build clusters of loci from reference and alternative annotations.
    
    Args:
        read_ref: Dictionary mapping chromosome_strand to list of Locus objects from reference annotation
        read_alt: Dictionary mapping chromosome_strand to list of Locus objects from alternative annotation
        dna_mol: Chromosome_strand key to process
        
    Returns:
        List of Cluster objects
    """
    cluster_list = []
    
    # Get the max end positions of the loci from both annotations
    sentinel = 1 + max(
        max(locus.end for locus in ref_loci) if ref_loci else 0,
        max(locus.end for locus in alt_loci) if alt_loci else 0
    )

    ref_sentinel = Locus.sentinel("sentinel_ref", sentinel)
    alt_sentinel = Locus.sentinel("sentinel_alt", sentinel)

    ref_loci.append(ref_sentinel)
    alt_loci.append(alt_sentinel)
    
    ref_i, alt_i = 0, 0
    ref_locus = ref_loci[ref_i]
    alt_locus = alt_loci[alt_i]
    
    cluster_max = -1
    cluster_loci = {"ref": [], "alt": []}
    cluster_id = 0
    
    while ref_locus.start < sentinel or alt_locus.start < sentinel:
        if ref_locus.start <= alt_locus.start:
            locus = ref_locus
            type = "ref"
            ref_locus, ref_i = ref_loci[ref_i+1], ref_i + 1
        else:
            locus = alt_locus
            type = "alt"
            alt_locus, alt_i = alt_loci[alt_i+1], alt_i + 1
        
        if locus.start > cluster_max:
            if cluster_max != -1:
                cluster_name = f"{dna_mol}_{cluster_id}"
                cluster_list.append(Cluster(cluster_name, cluster_loci["ref"], cluster_loci["alt"], cluster_max))
            cluster_id += 1
            cluster_loci = {"ref": [], "alt": []}
            cluster_max = -1
        
        cluster_max = max(cluster_max, locus.end)
        cluster_loci[type].append(locus)
    
    if cluster_max != -1:
        cluster_name = f"{dna_mol}_{cluster_id}"
        cluster_list.append(Cluster(cluster_name, cluster_loci["ref"], cluster_loci["alt"]  , cluster_max))

    ref_loci.clear()
    alt_loci.clear()
    return cluster_list

        