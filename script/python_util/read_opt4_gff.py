from array import array
from typing import Optional
from collections import defaultdict
import re
from locus import Locus

STRING_CACHE_DIRECT = "direct"
STRING_CACHE_REVERSE = "reverse"


def compile_id_regex():
    return re.compile(r'\bID=([^;\n\r]+)')

def compile_parent_regex():
    return re.compile(r'\bParent=([^;\n\r]+)')

def gff_to_cdsInfo(gff_file: str, relevant_gene_ids: Optional[set[str]] = None) -> dict[str, Locus]:
    """
    Optimized GFF parser: uses flat arrays and defaultdicts to reduce memory and improve speed.

    Args:
        gff_file: Path to the GFF file
        relevant_gene_ids: Optional set of gene IDs to filter (if None, all genes are processed)

    Returns:
        A dictionary mapping chromosome_strand keys to lists of Locus objects
    """
    id_regex = compile_id_regex()
    parent_regex = compile_parent_regex()

    mrna_id2CDS = defaultdict(lambda: array('L'))  # dict[str, array('L')]
    gene_id2mRNA = defaultdict(list)  # dict[str, list[str]]

    gene_ids = []
    gene_chr_ids = []
    gene_starts = array("L")
    gene_ends = array("L")
    gene_strands = array("B")  # 1 for '+', 0 for '-'

    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue
            infos = line.strip().split("\t")
            type = infos[2]

            if type == "CDS":
                start, end = int(infos[3]), int(infos[4])
                mrna_id = parent_regex.search(infos[8]).group(1)
                phase = int(infos[7]) if infos[7] != "." else 0
                mrna_id2CDS[mrna_id].extend([phase, start, end])

            elif type == "mRNA":
                mrna_id = id_regex.search(infos[8]).group(1)
                gene_id = parent_regex.search(infos[8]).group(1)

                if relevant_gene_ids and gene_id not in relevant_gene_ids:
                    continue

                gene_id2mRNA[gene_id].append(mrna_id)

            elif type == "gene":
                gene_id = id_regex.search(infos[8]).group(1)

                if relevant_gene_ids and gene_id not in relevant_gene_ids:
                    continue

                gene_ids.append(gene_id)
                gene_chr_ids.append(infos[0])
                gene_starts.append(int(infos[3]))
                gene_ends.append(int(infos[4]))
                gene_strands.append(1 if infos[6] == "+" else 0)

    chrStrand_2_loci = {}
    print(f"Number of genes (opt): {len(gene_ids)}")
    input("Press Enter to continue...")
    for i, gene_id in enumerate(gene_ids):
        chr_id = gene_chr_ids[i]
        start = gene_starts[i]
        end = gene_ends[i]
        strand_bool = gene_strands[i]
        strand = "+" if strand_bool else "-"

        gene_mrnas = gene_id2mRNA.pop(gene_id, [])
        if not gene_mrnas:
            continue

        mrna_dict = {}
        phases_dict = {}

        for mrna_id in gene_mrnas:
            if mrna_id in mrna_id2CDS:
                flat = mrna_id2CDS.pop(mrna_id)
                cds_list = [(flat[j], flat[j+1], flat[j+2]) for j in range(0, len(flat), 3)]
                cds_list.sort(key=lambda x: x[1])

                for j in range(len(cds_list) - 1):
                    if cds_list[j][2] >= cds_list[j+1][1]:
                        print(f"Error: mRNA {mrna_id} has overlapping CDS regions")
                        exit(1)

                mrna_dict[mrna_id] = array('L', (coord for cds in cds_list for coord in (cds[1], cds[2])))
                phases_dict[mrna_id] = cds_list[0][0] if strand_bool else cds_list[-1][0]

        if mrna_dict:
            locus = Locus(
                name=gene_id,
                start=start,
                end=end,
                mRNAs=mrna_dict,
                phases=phases_dict
            )
            direction = STRING_CACHE_DIRECT if strand_bool else STRING_CACHE_REVERSE
            chr_strand = f"{chr_id}_{direction}"
            chrStrand_2_loci.setdefault(chr_strand, []).append(locus)

    for chr_strand in chrStrand_2_loci:
        chrStrand_2_loci[chr_strand].sort(key=lambda l: (l.start, l.end))
    input("convertin genes done ")
    return chrStrand_2_loci
