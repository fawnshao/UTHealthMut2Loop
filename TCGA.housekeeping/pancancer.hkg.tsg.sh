# find hkg and tsg in TCGA pancancer data 
# use cohort.TCGA.Pan-Cancer.PANCAN log2(RSEM + 1)
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
input=tcga_RSEM_gene_tpm
# head -1 PANCAN_clinicalMatrix | cut -f 1,21-22,25 > $input.samples
# perl $myperl <(cut -f 1,21-22,25 PANCAN_clinicalMatrix) <(head_line $input | awk '{print $2}' | tail -n +2) 0 0 | cut -f 1,3- >> $input.samples
# awk -F"\t" '{print $1"\t"$3"|"$4}' $input.samples > $input.samples.sim
