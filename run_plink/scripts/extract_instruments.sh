#!/bin/bash

# Requires:
# bgenix
# plink2
# cat-bgen

outname="mrbase_instruments_20180612"

bgen_pattern="/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/ukb_imp_chrCHROM_v2.bgen"
bgen_index_pattern="/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/ukb_bgi_chrCHROM_v2.bgi"
snp_list="${HOME}/mr-eve/ukbb-instruments/run_plink/genotypes/${outname}.txt"

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

plink2 \
--bgen ../genotypes/${outname}.bgen \
--sample /mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/id_mapping/app8786/data.sample \
--hard-call-threshold 0.4999 \
--make-bed \
--out ../genotypes/${outname}

plink2 \
--bfile ../genotypes/${outname} \
--remove ../other_files/exclude.txt \
--freq \
--out ../genotypes/${outname}


# module load apps/qctool-2.0rc4

# qctool -g ${outname}.bgen -s ~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/8786_link/imp/ukb878_imp_chr22_v2_s487406.sample -excl-samples rem.txt -og instruments_fil.bgen 


