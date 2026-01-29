# -*- coding: utf-8 -*-

from pathlib import Path
from cdscompare.python_util.comparison import annotation_comparison
from cdscompare.python_util.annotation import Annotation, AnnotationPair

BASE = Path(__file__).parent.parent
DATA = BASE / ".." / "data" / "tests"

def make_pair(ref_path: Path, alt_path: Path) -> AnnotationPair:
   """Create an AnnotationPair from two GFF paths (test helper)."""
   ref = Annotation(id=ref_path.stem, path=ref_path)
   alt = Annotation(id=alt_path.stem, path=alt_path)
   return AnnotationPair(ref=ref, alt=alt)

# test function for the 'CDScompare.py' function 'annotation_comparison' ('main' function of the program)
def test_new_annotation_comparison(tmp_path: Path):

   # dictionary of inputs and expected ouputs for each test file for the 'CDScompare.py' function 'annotation_comparison' ('main' function of the program)
   test_dict = {
      "basic" : [DATA / "basic_test.gff3",
                  DATA / "identical_test.gff3",
                     {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                        "reference start" : 100,
                        "reference end" : 299,
                        "alternative" : 'chr2A_00611930',
                        "alternative start" : 100,
                        "alternative end" : 299,
                        "reference mRNA" : 'chr2A_00611930_mrna',
                        "alternative mRNA" : 'chr2A_00611930_mrna',
                        "mismatch/match" : [150, 0, 0],
                        "identity" : 100.0,
                        "mismatch zones" : ([], []),
                        "cluster name" : "chr2A_direct_1",
                        "reference mRNA number" : 1,
                        "alternative mRNA number" : 1}]]}],

      "minus-CDS" : [DATA / "basic_test.gff3",
                     DATA / "minus-CDS_test.gff3",
                     {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                        "reference start" : 100,
                        "reference end" : 299,
                        "alternative" : 'chr2A_00611930',
                        "alternative start" : 100,
                        "alternative end" : 299,
                        "reference mRNA" : 'chr2A_00611930_mrna',
                        "alternative mRNA" : 'chr2A_00611930_mrna',
                        "mismatch/match" : [90, 60, 0],
                        "identity" : 60.0,
                        "mismatch zones" : ([150, 209], []),
                        "cluster name" : "chr2A_direct_1",
                        "reference mRNA number" : 1,
                        "alternative mRNA number" : 1}]]}],

      "fusion" : [DATA / "basic_test.gff3",
                  DATA / "fusion_test.gff3",
                  {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                     "reference start" : 100,
                     "reference end" : 299,
                     "alternative" : 'chr2A_00611930',
                     "alternative start" : 100,
                     "alternative end" : 299,
                     "reference mRNA" : 'chr2A_00611930_mrna',
                     "alternative mRNA" : 'chr2A_00611930_mrna',
                     "mismatch/match" : [30, 20, 120],
                     "identity" : 17.65,
                     "mismatch zones" : ([130, 149], [150, 209, 240, 299]),
                     "cluster name" : "chr2A_direct_1",
                     "reference mRNA number" : 1,
                     "alternative mRNA number" : 1}]]}],

      "shift" : [DATA / "basic_test.gff3",
                  DATA / "shift_test.gff3",
                  {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                     "reference start" : 100,
                     "reference end" : 299,
                     "alternative" : 'chr2A_00611930',
                     "alternative start" : 100,
                     "alternative end" : 299,
                     "reference mRNA" : 'chr2A_00611930_mrna',
                     "alternative mRNA" : 'chr2A_00611930_mrna',
                     "mismatch/match" : [30, 1, 119],
                     "identity" : 20.0,
                     "mismatch zones" : ([150, 150], [151, 209, 240, 299]),
                     "cluster name" : "chr2A_direct_1",
                     "reference mRNA number" : 1,
                     "alternative mRNA number" : 1}]]}],

      "reverse-reverse" : [DATA / "reverse_test.gff3",
                  DATA / "reverse_test.gff3",
                  {'chr2A_reverse': [[{"reference" : 'chr2A_00611930',
                     "reference start" : 100,
                     "reference end" : 299,
                     "alternative" : 'chr2A_00611930',
                     "alternative start" : 100,
                     "alternative end" : 299,
                     "reference mRNA" : 'chr2A_00611930_mrna',
                     "alternative mRNA" : 'chr2A_00611930_mrna',
                     "mismatch/match" : [150, 0, 0],
                     "identity" : 100.0,
                     "mismatch zones" : ([], []),
                     "cluster name" : "chr2A_reverse_1",
                     "reference mRNA number" : 1,
                     "alternative mRNA number" : 1}]]}],

      "reverse-modified" : [DATA / "reverse_test.gff3",
                  DATA / "reverse_modified_test.gff3",
                  {'chr2A_reverse': [[{"reference" : 'chr2A_00611930',
                     "reference start" : 100,
                     "reference end" : 299,
                     "alternative" : 'chr2A_00611930',
                     "alternative start" : 100,
                     "alternative end" : 299,
                     "reference mRNA" : 'chr2A_00611930_mrna',
                     "alternative mRNA" : 'chr2A_00611930_mrna',
                     "mismatch/match" : [60, 10, 90],
                     "identity" : 37.5,
                     "mismatch zones" : ([210, 219], [100, 129, 150, 209]),
                     "cluster name" : "chr2A_reverse_1",
                     "reference mRNA number" : 1,
                     "alternative mRNA number" : 1}]]}],

      "diff-start-before" : [DATA / "basic_test.gff3",
                              DATA / "diff-start-before_test.gff3",
                              {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 40,
                                 "alternative end" : 239,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [30, 180, 30],
                                 "identity" : 12.5,
                                 "mismatch zones" : ([40, 69, 90, 99, 130, 179, 210, 299], [100, 129]),
                                 "cluster name" : "chr2A_direct_1",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}]]}],

      "diff-start-after" : [DATA / "basic_test.gff3",
                           DATA / "diff-start-after_test.gff3",
                           {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                              "reference start" : 100,
                              "reference end" : 299,
                              "alternative" : 'chr2A_00611930',
                              "alternative start" : 160,
                              "alternative end" : 359,
                              "reference mRNA" : 'chr2A_00611930_mrna',
                              "alternative mRNA" : 'chr2A_00611930_mrna',
                              "mismatch/match" : [30, 180, 30],
                              "identity" : 12.5,
                              "mismatch zones" : ([100, 129, 150, 159, 190, 239, 270, 359], [160, 189]),
                              "cluster name" : "chr2A_direct_1",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1}]]}],

      "basic-2-loci" : [DATA / "basic-2-loci_test.gff3",
                        DATA / "identical-2-loci_test.gff3",
                           {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                              "reference start" : 100,
                              "reference end" : 299,
                              "alternative" : 'chr2A_00611930',
                              "alternative start" : 100,
                              "alternative end" : 299,
                              "reference mRNA" : 'chr2A_00611930_mrna',
                              "alternative mRNA" : 'chr2A_00611930_mrna',
                              "mismatch/match" : [150, 0, 0],
                              "identity" : 100.0,
                              "mismatch zones" : ([], []),
                              "cluster name" : "chr2A_direct_1",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1}],
                              [{"reference" : 'chr2A_00620000',
                                 "reference start" : 600,
                                 "reference end" : 899,
                                 "alternative" : 'chr2A_00620000',
                                 "alternative start" : 600,
                                 "alternative end" : 899,
                                 "reference mRNA" : 'chr2A_00620000_mrna',
                                 "alternative mRNA" : 'chr2A_00620000_mrna',
                                 "mismatch/match" : [200, 0, 0],
                                 "identity" : 100.0,
                                 "mismatch zones" : ([], []),
                                 "cluster name" : "chr2A_direct_2",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}]]}],

      "minus-loci" : [DATA / "basic_test.gff3",
                        DATA / "identical-2-loci_test.gff3",
                           {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                              "reference start" : 100,
                              "reference end" : 299,
                              "alternative" : 'chr2A_00611930',
                              "alternative start" : 100,
                              "alternative end" : 299,
                              "reference mRNA" : 'chr2A_00611930_mrna',
                              "alternative mRNA" : 'chr2A_00611930_mrna',
                              "mismatch/match" : [150, 0, 0],
                              "identity" : 100.0,
                              "mismatch zones" : ([], []),
                              "cluster name" : "chr2A_direct_1",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1}],
                              [{"reference" : '~',
                                 "reference start" : '_',
                                 "reference end" : '_',
                                 "alternative" : 'chr2A_00620000',
                                 "alternative start" : 600,
                                 "alternative end" : 899,
                                 "reference mRNA" : '_',
                                 "alternative mRNA" : 'chr2A_00620000_mrna',
                                 "mismatch/match" : [],
                                 "identity" : 0.0,
                                 "mismatch zones" : '_',
                                 "cluster name" : "chr2A_direct_2",
                                 "reference mRNA number" : '_',
                                 "alternative mRNA number" : 1}]]}],

      "overlapping-loci" : [DATA / "overlapping-loci_test.gff3",
                           DATA / "overlapping-loci-alt_test.gff3",
                           {'chr2A_direct': [[{'reference': 'chr2A_1000',
                              'reference start': 50,
                              'reference end': 149,
                              'alternative': 'chr2A_1000',
                              'alternative start': 100,
                              'alternative end': 249,
                              "reference mRNA" : 'chr2A_1000_mrna',
                              "alternative mRNA" : 'chr2A_1000_mrna',
                              'mismatch/match': [0, 150, 50],
                              'identity': 0.0,
                              'mismatch zones': ([50, 99, 150, 249], [100, 149]),
                              "cluster name" : "chr2A_direct_1",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1},
                              {'reference': 'chr2A_2000',
                              'reference start': 200,
                              'reference end': 349,
                              'alternative': 'chr2A_2000',
                              'alternative start': 300,
                              'alternative end': 449,
                              "reference mRNA" : 'chr2A_2000_mrna',
                              "alternative mRNA" : 'chr2A_2000_mrna',
                              'mismatch/match': [0, 200, 50],
                              'identity': 0.0,
                              'mismatch zones': ([200, 299, 350, 449], [300, 349]),
                              "cluster name" : "chr2A_direct_1",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1},
                              {'reference': 'chr2A_3000',
                              'reference start': 400,
                              'reference end': 549,
                              'alternative': 'chr2A_3000',
                              'alternative start': 500,
                              'alternative end': 599,
                              "reference mRNA" : 'chr2A_3000_mrna',
                              "alternative mRNA" : 'chr2A_3000_mrna',
                              'mismatch/match': [0, 150, 50],
                              'identity': 0.0,
                              'mismatch zones': ([400, 499, 550, 599], [500, 549]),
                              "cluster name" : "chr2A_direct_1",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1}],
                              [{'reference': 'chr2A_4000',
                              'reference start': 650,
                              'reference end': 699,
                              'alternative': 'chr2A_4000',
                              'alternative start': 650,
                              'alternative end': 699,
                              "reference mRNA" : 'chr2A_4000_mrna',
                              "alternative mRNA" : 'chr2A_4000_mrna',
                              'mismatch/match': [50, 0, 0],
                              'identity': 100.0,
                              'mismatch zones': ([], []),
                              "cluster name" : "chr2A_direct_2",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1}],
                              [{'reference': 'chr2A_5000',
                              'reference start': 750,
                              'reference end': 799,
                              'alternative': 'chr2A_5000',
                              'alternative start': 750,
                              'alternative end': 779,
                              "reference mRNA" : 'chr2A_5000_mrna',
                              "alternative mRNA" : 'chr2A_5000_mrna',
                              'mismatch/match': [30, 20, 0],
                              'identity': 60.0,
                              'mismatch zones': ([780, 799], []),
                              "cluster name" : "chr2A_direct_3",
                              "reference mRNA number" : 1,
                              "alternative mRNA number" : 1},
                              {'reference': '~',
                              'reference start': '_',
                              'reference end': '_',
                              'alternative': 'chr2A_6000',
                              'alternative start': 790,
                              'alternative end': 849,
                              "reference mRNA" : '_',
                              "alternative mRNA" : 'chr2A_6000_mrna', # if no pairing with reference no mRNA picking
                              #"alternative mRNA" : '-',
                              'mismatch/match': [],
                              'identity': 0.0,
                              'mismatch zones': '_',
                              "cluster name" : "chr2A_direct_3",
                              "reference mRNA number" : '_',
                              "alternative mRNA number" : 1}]]}],

      'length_computation' : [DATA / "length_computation_ref_test.gff3",
                              DATA / "length_computation_alt_test.gff3",
                              {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 8,
                                 "reference end" : 13,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 1,
                                 "alternative end" : 13,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [2, 8, 0],
                                 "identity" : 20.0,
                                 "mismatch zones" : ([1, 4, 8, 11], []),
                                 "cluster name" : "chr2A_direct_1",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}]]}],
   
      "multi_chromosome" : [DATA / "multi_chromosome_test_ref.gff3",
                  DATA / "multi_chromosome_test_alt.gff3",
                  {'chr2A_direct': [[{"reference" : 'chr2A_00611930',
                     "reference start" : 100,
                     "reference end" : 299,
                     "alternative" : 'chr2A_00611930',
                     "alternative start" : 100,
                     "alternative end" : 299,
                     "reference mRNA" : 'chr2A_00611930_mrna',
                     "alternative mRNA" : 'chr2A_00611930_mrna',
                     "mismatch/match" : [30, 1, 119],
                     "identity" : 20.0,
                     "mismatch zones" : ([150, 150], [151, 209, 240, 299]),
                     "cluster name" : "chr2A_direct_1",
                     "reference mRNA number" : 1,
                     "alternative mRNA number" : 1}]],
                  'chr2B_direct': [[{"reference" : 'chr2B_00620000',
                     "reference start" : 600,
                     "reference end" : 899,
                     "alternative" : 'chr2B_00620000',
                     "alternative start" : 600,
                     "alternative end" : 899,
                     "reference mRNA" : 'chr2B_00620000_mrna',
                     "alternative mRNA" : 'chr2B_00620000_mrna',
                     "mismatch/match" : [200, 0, 0],
                     "identity" : 100.0,
                     "mismatch zones" : ([], []),
                     "cluster name" : "chr2B_direct_1",
                     "reference mRNA number" : 1,
                     "alternative mRNA number" : 1}]]}],

      "diff_dnaMol1" : [DATA / "diff_dnaMol-1-chr_test.gff3",
                        DATA / "diff_dnaMol-2-chr_test.gff3",
                        {'chr2A_direct': [[{'reference': 'chr2A_00611930',
                              'reference start': 100,
                              'reference end': 299,
                              'alternative': 'chr2A_00611930',
                              'alternative start': 100,
                              'alternative end': 299,
                              'mismatch/match': [150, 0, 0],
                              'identity': 100.0,
                              'mismatch zones': ([], []),
                              'cluster name': 'chr2A_direct_1',
                              'reference mRNA': 'chr2A_00611930_mrna',
                              'alternative mRNA': 'chr2A_00611930_mrna',
                              'reference mRNA number': 1,
                              'alternative mRNA number': 1}]],
                        'chr2B_direct': [[{'reference': '~',
                              'reference start': '_',
                              'reference end': '_',
                              'alternative': 'chr2B_00620000',
                              'alternative start': 600,
                              'alternative end': 899,
                              'mismatch/match': [],
                              'identity': 0.0,
                              'mismatch zones': '_',
                              'cluster name': 'chr2B_direct_1',
                              'reference mRNA': '_',
                              'alternative mRNA': 'chr2B_00620000_mrna',
                              'reference mRNA number': '_',
                              'alternative mRNA number': 1}]]}],

      "diff_dnaMol2" : [DATA / "diff_dnaMol-2-chr_test.gff3",
                        DATA / "diff_dnaMol-1-chr_test.gff3",
                        {'chr2A_direct': [[{'reference': 'chr2A_00611930',
                              'reference start': 100,
                              'reference end': 299,
                              'alternative': 'chr2A_00611930',
                              'alternative start': 100,
                              'alternative end': 299,
                              'mismatch/match': [150, 0, 0],
                              'identity': 100.0,
                              'mismatch zones': ([], []),
                              'cluster name': 'chr2A_direct_1',
                              'reference mRNA': 'chr2A_00611930_mrna',
                              'alternative mRNA': 'chr2A_00611930_mrna',
                              'reference mRNA number': 1, 'alternative mRNA number': 1}]],
                        'chr2B_direct': [[{'reference': 'chr2B_00620000',
                              'reference start': 600,
                              'reference end': 899,
                              'alternative': '~',
                              'alternative start': '_',
                              'alternative end': '_',
                              'mismatch/match': [],
                              'identity': 0.0,
                              'mismatch zones': '_',
                              'cluster name': 'chr2B_direct_1',
                              'reference mRNA': 'chr2B_00620000_mrna',
                              'alternative mRNA': '_',
                              'reference mRNA number': 1,
                              'alternative mRNA number': '_'}]]}]


   }
   
   for name, (ref_path, alt_path, expected) in test_dict.items():
      out_dir = tmp_path / name
      out_dir.mkdir()
      pair = make_pair(ref_path, alt_path)
      result = annotation_comparison(pair, out_dir, True)
      assert result == expected, f"Mismatch in test case: {name}"

def test_geneoverlap_mRNA_dont(tmp_path: Path):
   gff_ref=DATA / "cluster_3274_urgi.gff"
   gff_alt=DATA / "cluster_3274_ncbi.gff"
   pair = make_pair(gff_ref, gff_alt)
   annotation_comparison(pair, tmp_path, True)

def test_missed_genomic_overlap(tmp_path: Path):   
   gff_ref=DATA / "missed_overlap_ref.gff"
   gff_alt=DATA / "missed_overlap_alt.gff"
   pair = make_pair(gff_ref, gff_alt)
   annotation_comparison(pair, tmp_path, True)

def test_phase(tmp_path: Path):   
   gff_ref=DATA / "phased2_ref.gff"
   gff_alt=DATA / "phased2_alt.gff"
   pair = make_pair(gff_ref, gff_alt)
   result = annotation_comparison(pair, tmp_path, True)
   match = result["Chr6D_direct"][0][0]
   identity = match['identity']
   assert (identity>90)

def test_phase_rev(tmp_path: Path):   
   gff_ref=DATA / "phased-rev_ref.gff"
   gff_alt=DATA / "phased-rev_alt.gff"
   pair = make_pair(gff_ref, gff_alt)
   result = annotation_comparison(pair, tmp_path, True)
   match = result["Chr6D_reverse"][0][0]
   identity = match['identity']
   assert (identity>60)

# def test_multi_mrnas():   
#    with tempfile.TemporaryDirectory() as tmpdir:
#       gff_ref="data/tests/multi_mrna_BW_ref.gff"
#       gff_alt="data/tests/multi_mrna_BW_alt.gff"

#    print(tmpdir)
#    result = annotation_comparison(gff_ref, gff_alt, tmpdir, True)
#    match = result["Chr7D_reverse"][0][0]
#    identity = match['identity']
#    assert (identity>90)

# def test_full():   
#    with tempfile.TemporaryDirectory() as tmpdir:
#       gff_ref="/Users/ranwez/My_data/Projects/2023_ANNOT_WHEAT_LRR/BleTendre/01_sorted_input_gffs/ref_chr_urgi_HC_unsorted.gff"
#       gff_alt="/Users/ranwez/My_data/Projects/2023_ANNOT_WHEAT_LRR/BleTendre/01_sorted_input_gffs/alt_chr_ncbi_sorted.gff"
      

#    print(tmpdir)
#    result = annotation_comparison(gff_ref, gff_alt, tmpdir, True)
#    match = result["Chr6D_direct"][0][0]
#    identity = match['identity']
#    assert (identity>90)
