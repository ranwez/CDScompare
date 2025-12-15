from array import array
from typing import Optional
from collections import defaultdict
import python_util.intervals as iu
import re
 

STRING_CACHE_DIRECT = "direct"
STRING_CACHE_REVERSE = "reverse"

parent_regex = re.compile(r'Parent=([^;\s]+)')
id_regex = re.compile(r'ID=([^;\s]+)')
id_parent_regex = re.compile(r'ID=([^;\s]+).*Parent=([^;\s]+)')

def extract_id(attr_str: str) -> Optional[str]:
    for part in attr_str.split(';'):
        part = part.strip()
        if part.startswith("ID="):
            return part[3:].strip()
    return None

def extract_parent_old(attr_str: str) -> Optional[str]:
    for part in attr_str.split(';'):
        part = part.strip()
        if part.startswith("Parent="):
            return part[7:].strip()
    return None
def extract_parent(attr_str: str) -> Optional[str]:
    return parent_regex.search(attr_str).group(1) 


def extract_parent_id(attr_str: str) -> tuple[Optional[str], Optional[str]]:
    id_val = parent_val = None
    for part in attr_str.split(';'):
        part = part.strip()
        if part.startswith("ID="):
            id_val = part[3:].strip()
            if parent_val:
                break
        elif part.startswith("Parent="):
            parent_val = part[7:].strip()
            if id_val :
                break
    return id_val, parent_val

class Locus:
    """
    Represents an annotation's locus identified from a GFF file.
    Uses __slots__ to reduce memory footprint.
    """
    __slots__ = ('name', 'mRNAs', 'start', 'end', 'phases', 'ids', 'mrna_intervals')

    def __init__(self, name="", mRNAs=None, phases=None, ids=None, start=-1, end=-1):
        """
        Initialize a Locus object.

        Args:
            name: ID of the locus
            mRNAs: List of CDS coordinate arrays (array('L'))
            phases: array('B') of CDS phases
            ids: list of mRNA IDs
            start: Start coordinate of the locus
            end: End coordinate of the locus
        """
        self.name = name
        self.mRNAs = mRNAs if mRNAs is not None else []
        self.phases = phases if phases is not None else array('B')
        self.ids = ids if ids is not None else []
        self.start = start
        self.end = end
        self.mrna_intervals = None
    def set_mrna_intervals(self):
        self.mrna_intervals = [iu.OrderedIntervals(mRNA, True) for mRNA in self.mRNAs]

    @classmethod
    def sentinel(cls, name: str, end: int) -> "Locus":
        return cls(
            name=name,
            mRNAs=[],
            phases=array("B"),
            ids=[],
            start=end,
            end=end + 1
        )
        
    @classmethod
    def builder(cls, name: str = "", mRNAs: dict = None, start: int = -1, end: int = -1, direction: str = "") -> "Locus":
        """
        Factory method to build a Locus object from dictionary-based mRNAs format.
        This is for backward compatibility with tests written for the older Locus implementation.
        
        Args:
            name: ID of the locus
            mRNAs: Dictionary mapping mRNA ID to a list of CDS coordinates
            start: Start coordinate of the locus
            end: End coordinate of the locus
            direction: Strand direction ('direct' or 'reverse')
            
        Returns:
            A new Locus instance
        """
        if mRNAs is None:
            mRNAs = {}
            
        # Convert dictionary to the new format
        mrna_list = []
        phase_list = array('B')
        mrna_ids = []
        
        for mrna_id, coords in mRNAs.items():
            mrna_list.append(array('L', coords))
            phase_list.append(0)  # Default phase
            mrna_ids.append(mrna_id)
            
        return cls(
            name=name,
            mRNAs=mrna_list,
            phases=phase_list,
            ids=tuple(mrna_ids),
            start=start,
            end=end
        )

    def reverse(self, cluster_end):
        """
        Reverse the coordinates of all mRNAs relative to cluster_end.
        Used for loci on the reverse strand.

        Args:
            cluster_end: End position of the parent cluster
        """
        for i in range(len(self.mRNAs)):
            mrna = self.mRNAs[i]
            length = len(mrna)
            
            for j in range(length // 2):
                front_val = cluster_end - mrna[length - j - 1]
                mrna[length - j - 1]  = cluster_end - mrna[j]
                mrna[j] = front_val
                

    def __str__(self):
        return f"Locus '{self.name}' ({self.start}-{self.end})"

    def show_init(self):
        return f"Locus(name='{self.name}', mRNAs={self.mRNAs}, start={self.start}, end={self.end}')"

def gff_to_cdsInfo(gff_file: str) -> dict[str, list[Locus]]:
    """
    Optimized GFF parser: uses flat arrays and defaultdicts to reduce memory and improve speed.

    Args:
        gff_file: Path to the GFF file
        relevant_gene_ids: Optional set of gene IDs to filter (if None, all genes are processed)

    Returns:
        A dictionary mapping chromosome_strand keys to lists of Locus objects
    """

    mrna_id2CDS = defaultdict(lambda: array('L'))  # dict[str, array('L')]
    gene_id2mRNA = defaultdict(list)  # dict[str, list[str]]

    gene_ids = []
    gene_chr_ids = []
    gene_starts = array("L")
    gene_ends = array("L")
    gene_strands = array("B")  # 1 for '+', 0 for '-'

    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#") or line.isspace():
                continue
            infos = line.split("\t")
            type = infos[2]
            if type == "CDS":
                start, end = int(infos[3]), int(infos[4])
                mrna_id = parent_regex.search(infos[8]).group(1)
                phase = int(infos[7]) if infos[7] != "." else 0
                mrna_id2CDS[mrna_id].extend([phase, start, end])

            elif type == "mRNA":
                id_parent = id_parent_regex.search(infos[8])
                mrna_id = id_parent.group(1)
                gene_id = id_parent.group(2)
                gene_id2mRNA[gene_id].append(mrna_id)

            elif type == "gene":
                gene_id = id_regex.search(infos[8]).group(1)
                gene_ids.append(gene_id)
                gene_chr_ids.append(infos[0])
                gene_starts.append(int(infos[3]))
                gene_ends.append(int(infos[4]))
                gene_strands.append(1 if infos[6] == "+" else 0)
    chrStrand_2_loci = defaultdict(list) 
    print(f"Converting {len(gene_ids)} genes to loci...")
    #input("Press Enter to continue...")
    empty=[]
    for i, gene_id in enumerate(gene_ids):

        gene_mrnas = gene_id2mRNA.pop(gene_id, empty)
        if not gene_mrnas:
            continue
        chr_id = gene_chr_ids[i]
        start = gene_starts[i]
        end = gene_ends[i]
        strand_bool = gene_strands[i]

        mrna_list = []
        phase_list = array('B')
        mrna_ids = []
        for mrna_id in gene_mrnas:
            if mrna_id in mrna_id2CDS:
                flat = mrna_id2CDS.pop(mrna_id)
                cds_list = sorted(
                   ((flat[j], flat[j+1], flat[j+2]) for j in range(0, len(flat), 3)),
                    key=lambda x: x[1]
                )

                for j in range(len(cds_list) - 1):
                    if cds_list[j][2] >= cds_list[j+1][1]:
                        print(f"Error: mRNA {mrna_id} has overlapping CDS regions")
                        exit(1)

                mrna_list.append(array('L', (coord for cds in cds_list for coord in (cds[1], cds[2]))))
                phase_list.append(cds_list[0][0] if strand_bool else cds_list[-1][0])
                mrna_ids.append(mrna_id)

        if mrna_list:
            locus = Locus(
                name=gene_id,
                mRNAs=mrna_list,
                phases=phase_list,
                ids=tuple(mrna_ids),
                start=start,
                end=end
            )
            direction = STRING_CACHE_DIRECT if strand_bool else STRING_CACHE_REVERSE
            chr_strand = f"{chr_id}_{direction}"
            #chrStrand_2_loci.setdefault(chr_strand, []).append(locus)
            chrStrand_2_loci[chr_strand].append(locus)
    for chr_strand in chrStrand_2_loci:
        chrStrand_2_loci[chr_strand].sort(key=lambda l: (l.start, l.end))

    return chrStrand_2_loci

