import pandas as pd

import argparse

# "/Users/morril01/Documents/PhD/CDA_in_Cancer/out/ProjectSubtractingExposures/out_ffad9288-c622-11e3-bf01-24c6515278c0.consensus.20160830.somatic.snv_mnv.vcf"
# "/Users/morril01/Documents/PhD/CDA_in_Cancer/data/pcawg/consensus_subclonal_reconstruction_mutccf_20170325/ffad9288-c622-11e3-bf01-24c6515278c0_mutation_ccf.txt"

parser = argparse.ArgumentParser(description='Process input file')
parser.add_argument('--file_VAF', type=str,
                    help='file to process (1/2)')
parser.add_argument('--file_ccf', type=str,
                    help='file to process (2/2)')
parser.add_argument('--out_file', type=str,
                    help='output file')
args = parser.parse_args()


DF1=pd.read_csv(args.file_VAF, names=['chromosome', 'position', 'flanking1', 'flanking2', 'mutation', 'VAF'],
	comment='#', delim_whitespace=True)
DF2=pd.read_csv(args.file_ccf, header=0, sep='\t')

# print(DF1)
# print(DF2)

merged_df = DF2.merge(DF1, how = 'inner', on = ['chromosome', 'position'])
# print(merged_df)
merged_df.to_csv(args.out_file, sep="\t")