#!/bin/bash

# Requires:
# bgenix
# plink2
# cat-bgen

bgen_pattern="${HOME}/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_imp_chrCHROM_v2.bgen"
bgen_index_pattern="${HOME}/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_bgi_chrCHROM_v2.bgi"
snp_list="${HOME}/repo/ukbb-instruments/run_plink/genotypes/mrbase_instruments_20180612.txt"
outname="mrbase_instruments_20180612"

# Extract using bgen
temp_prefix="../genotypes/temp_genotypes"
for chrom in {1..22}; do
    inbgen=${bgen_pattern/CHROM/$chrom}
    inbgenidx=${bgen_index_pattern/CHROM/$chrom}
    bgenix -g ${inbgen} -i ${inbgenidx} -incl-rsids ${snp_list} > ${temp_prefix}.${chrom}.bgen
done

cmd=""
for chrom in {1..22}; do
	cmd="${cmd} ${temp_prefix}.${chrom}.bgen"
done

cat-bgen -g ${cmd} -og ../genotypes/${outname}.bgen
cp ~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/id_mapping/app8786/data.sample ../genotypes/app8786.sample

plink2 \
--bgen ../genotypes/${outname}.bgen \
--sample ~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/id_mapping/app8786/data.sample \
--hard-call-threshold 0.4999 \
--make-bed \
--out ../genotypes/${outname}

plink2 \
--bfile ../genotypes/${outname} \
--freq \
--out ../genotypes/${outname}

# module load apps/qctool-2.0rc4

# qctool -g instruments.bgen -s ~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/8786_link/imp/ukb878_imp_chr22_v2_s487406.sample -excl-samples rem.txt -og instruments_fil.bgen 


