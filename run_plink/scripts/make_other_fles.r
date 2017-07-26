# For bgenie

# exclusion lists - recommended + ancestry (stringent)
# sample file / fam file
# read in file, map to 15825, subset and reorder
# Set excluded to NA
# Set NA to -999
convert_mapping <- function(idlist, mapping)
{
	index <- match(idlist, mapping$app15825)
	out <- mapping$app8786[index]
	return(out)
}

reorder_table <- function(x, samp, idname)
{
	a <- left_join(samp, x, by=c("id"=idname))
	a <- a[order(a$index), ]
	a <- select(a, -index)
	return(a)
}


library(data.table)
library(dplyr)
path <- "~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/derived/"

# Get mappings
mapping <- fread("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/data.7445.csv", he=T, sep=",")

# Get exclusion list
rel <- paste0(path, "related_ids/data.relateds_exclusions.txt")
rel <- as.numeric(fread(rel)[[1]])
exc <- paste0(path, "standard_exclusions/meta.recommended_exclusions.txt")
exc <- as.numeric(fread(exc)[[1]])
anc <- paste0(path, "ancestry/data.non_europeans_exclusions.txt")
anc <- as.numeric(fread(anc)[[1]])

to_remove <- unique(c(rel, exc, anc))


# Read in the sample file
samp <- "~/data/bb/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/8786_link/imp/ukb878_imp_chr22_v2_s487406.excl-samples"
samp <- fread(samp, skip=2)
samp <- data.frame(id=samp$V1, index=1:nrow(samp))


# Make covariate file
covs <- paste0(path, "standard_covariates/data.covariates.txt")
covs <- fread(covs)
covs <- covs[,c("V2"):=NULL]
covs$V3[covs$V3 == "M"] <- 1
covs$V3[covs$V3 == "F"] <- 2
table(covs$V3)
covs$V4[covs$V4 == "UKBB"] <- 1
covs$V4[covs$V4 == "UKBL"] <- 2

covs2 <- reorder_table(covs, samp, "V1")

pcs <- paste0(path, "principal_components/data.pca1-10.txt")
pcs <- fread(pcs)
pcs <- pcs[, c("V2"):=NULL]
covs <- left_join(covs2, pcs, by=c("id"="V1"))

stopifnot(all(covs$id == samp$id))
covs$id <- paste(covs$id, covs$id)

write.table(covs, file="../other_files/covs.txt", row=F, col=F, qu=F)
write.table(data.frame(to_remove, to_remove), file="../other_files/exclude.txt", row=F, col=F, qu=F)
