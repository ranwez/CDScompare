
Description for each test file used in the program's unit tests:
(For a graphical overview, check the test_files.odp file) 

- basic_test : GFF file describing a unique locus with 3 CDS, intended to be used as a 'reference' for the other one-locus files

- identical_test : GFF file describing an identical locus to the one of basic_test, intended to test wether comparing an annotation file against itself yields 100% identity

- minus-CDS_test : GFF file describing a unique locus identical to basic_test except the second CDS is missing, which should yield an identity of 60% (90 matchs, 60 mismatchs).

- shift_test : GFF file describing the same locus as basic_test, except the start of the second codon is one nucleotide late, which should disrupt the codon position of each CDS nucleotide that follows, yielding an identity of 20% (30 matchs, 120 mismatches). It is intended to test wether the codon disruption follows on the next CDS and affects the final identity.

- fusion_test : GFF file describing the same locus as basic_test, except the first and second CDS have been fused (deleting the first intron of the locus). This disrupts the codon position of each CDS nucleotide that follows, yielding an identity of 17,6% (30 matchs, 140 mismatches). Since the fusion of structures is a type of error frequently observed in LRR annotation predictions, this file is intended to allow testing wether it affects the identity level.

- diff-start-before_test : GFF file describing the same locus as basic_test, except it begins at a position 60 nucleotides earlier than basic_test, yielding an identity of 12,5% (30 matches, 210 mismatches)

- diff-start-before_test : GFF file describing the same locus as basic_test, except it begins at a position 60 nucleotides later than basic_test, yielding an identity of 12,5% (30 matches, 210 mismatches)


- reverse_test : GFF file describing the same locus as basic_test but reversed, intended to test wether the program does assign an identity of 0% when comparing loci on different strands, and if the program handles well processing loci on the reverse strand.

- reverse_modified_test.gff3 : GFF file describing the same locus as reverse_test.gff3, but which presents a reading phase shift. It is intended to test if the program works well on reverse-strand loci. It should yield 37,5% identity (60 matches, 10+90 mismatches) when compared with reverse_test.gff3.


- multiple_mRNAs_test.gff3 : GFF file describing a single locus with multiple mRNAs. It is intended to test the selection of the maximum-identity mRNA from the locus' mRNAs. It should yield 80% identity (120 matches, 30 mismatches) when compared with basic_test.gff3.


- basic-2-loci_test : GFF file describing 2 loci with 3 and 2 CDS, intended to be used as a 'reference' for the other two-loci files. It an also be compared against basic_test to test wether the program does assign 100% to the first locus and 0% to the second one (absent in one but present in the other).

- identical-2-loci_test : GFF file describing the same 2 loci as basic-2-loci_test, intended to test wether the program does yield an identity of 100% when comparing two annotations with multiple loci.


- overlapping-loci_test : GFF file describing 5 one-CDS loci, intended to be used as a 'reference' when testing the program's gestion of overlapping reference-alternative loci 

- overlapping-loci-2_test : GFF file describing 6 loci, which when compared with overlapping-loci_test, should trigger the fusion of the first 3 loci of each annotation (as they are overlapping each other). The comparison of this file with overlapping-test_loci should yield the following identities: 0, 0, 0, 100, 60.


- length_computation_ref_test.gff3 : GFF file describing one mono-exonic loci, intended to be used as a reference when testing the program's comparison zone length computation.

- length_computation_alt_test.gff3 : GFF file describing one mono-exonic loci, which should yield an identity of 20% (8 mismatches, 2 matches) when compared with length_computation_ref_test.gff3.


- multi_chromosome_test_ref.gff3 : GFF file describing two loci on different chromosomes (chr2A and chr2B), intended to be used as reference when testing the program's chromosome gestion.

- multi_chromosome_test_alt.gff3 : GFF file describing two loci on different chromosomes (chr2A and chr2B), which should yield identities of 20 and 100 % when compared with multi_chromosome_test_ref.gff3.



