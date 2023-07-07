## /usr/bin/python3.6
## from https://www.biostars.org/p/334253/

import pysam
import argparse

parser = argparse.ArgumentParser(description='Process input file')
parser.add_argument('--file', type=str,
                    help='file to process')
args = parser.parse_args()

# open vcf file
vcf = pysam.VariantFile(args.file)
# open fasta file (ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz)
# genome = pysam.FastaFile("/Users/morril01/links_to_projects/subtractSignatures/other_files/genome.fa")
genome = pysam.FastaFile("/Users/morril01/Documents/PhD/CDA_in_Cancer/data/genome.fa")
# define by how many bases the variant should be flanked
flank = 1

# iterate over each variant
for record in vcf:
    # extract sequence
    #
    # The start position is calculated by subtract the number of bases
    # given by 'flank' from the variant position. The position in the vcf file
    # is 1-based. pysam's fetch() expected 0-base coordinate. That's why we
    # need to subtract on more base.
    #
    # The end position is calculated by adding the number of bases
    # given by 'flank' to the variant position. We also need to add the length
    # of the REF value and subtract again 1 due to the 0-based/1-based thing.
    #
    # Now we have the complete sequence like this:
    # [number of bases given by flank]+REF+[number of bases given by flank]

    ## the -1 was added on 20191209, because otherwise everything was one too separated (e.g. first nucleotide was 1, last trinucleotide was 4
    seq = genome.fetch(record.chrom, record.pos-1-flank, record.pos-1+len(record.ref)+flank)

    # print("{}\t{}\t{}".format(record.chrom, record.pos-1-flank, record.pos-1+len(record.ref)-1+flank))
    # print out tab seperated columns:
    # CRHOM, POS, REF, ALT, flanking sequencing with variant given in the format '[REF/ALT]'
    print(
        record.chrom,
        record.pos,
        record.ref,
        record.alts[0],
        '{}[{}/{}]{}'.format(seq[:flank], record.ref, record.alts[0], seq[flank+len(record.ref):]),
        sep="\t"
    )
