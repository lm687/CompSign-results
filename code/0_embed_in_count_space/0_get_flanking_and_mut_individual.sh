vcf_filename=$1
vcf_filename_tbi=$2 ## not actually used, but since the vcf needs to be indexed first, we use this as input for snakemake to infer the dependency
ccf_file=$3
out_file_final=$4

# echo 'File being processed is' $vcf_filename

# echo "Reading files"

raw_vcf_filename=`echo $vcf_filename | sed 's/.consensus.20160830.somatic.snv_mnv.vcf.gz//' | sed 's/pcawg_snv/consensus_subclonal_reconstruction_mutccf_20170325/'`

# echo "ccf file:\t" $ccf_file
# echo "raw_vcf_filename file:\t" $raw_vcf_filename

## get flanking bases
# echo "Running getflanking"
# echo ${out_file_final}tmpflanking
python3 0_embed_in_count_space/getflanking.py --file $vcf_filename > ${out_file_final}tmpflanking

# get the actual change
# echo "Running vcf-query"
# vcf-query -f '%POS,%QUAL,%INFO/VAF\n' $vcf_filename > ${out_file_final}tmpVAF
# cat ${out_file_final}tmpVAF |   sed 's/.*,//' > ${out_file_final}tmpVAF2
vcf-query -f '%POS,%QUAL,%INFO/VAF\n' $vcf_filename |  sed 's/.*,//' > ${out_file_final}tmpVAF2


## paste together
paste ${out_file_final}tmpflanking ${out_file_final}tmpVAF2 | column -s $'\t' -t > ${out_file_final}tmpout

# put ccf instead of vaf
# echo "Running match_ccf"
python3 0_embed_in_count_space/match_ccf.py --file_VAF ${out_file_final}tmpout --file_ccf ${ccf_file} --out_file $out_file_final

# echo 'File' $out_file_final 'created'

rm ${out_file_final}tmpflanking
rm ${out_file_final}tmpVAF2
rm  ${out_file_final}tmpout

