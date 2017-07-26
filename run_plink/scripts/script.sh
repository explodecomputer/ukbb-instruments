#!/bin/bash


bgen_pattern="${HOME}/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/ukb_imp_chrCHROM_v2.bgen"
bgen_index_pattern="${HOME}/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/ukb_bgi_chrCHROM_v2.bgi"
snp_list="${HOME}/repo/mr-eve/results/01/instrumentlist.txt"

# Extract using bgen
temp_prefix="temp_genotypes"
for chrom in {1..22}; do
    inbgen=${bgen_pattern/CHROM/$chrom}
    inbgenidx=${bgen_index_pattern/CHROM/$chrom}
    bgenix -g ${inbgen} -i ${inbgenidx} -incl-rsids ${snp_list} > ${temp_prefix}.${chrom}.bgen
done

cmd=""
for chrom in {1..22}; do
	cmd="${cmd} ${temp_prefix}.${chrom}.bgen"
done

cat-bgen -g ${cmd} -og instruments.bgen

plink2 \
--bgen instruments.bgen \
--sample ~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/8786_link/imp/ukb878_imp_chr22_v2_s487406.sample \
--hard-call-threshold 0.4999 \
--make-bed \
--out instrumentsc