#!/bin/bash

library(dplyr)
library(data.table)

plink_run <- function(geno, covs, pheno, pheno_name, out, exclude, fr, binary)
{
	flag <- ifelse(binary, " --logistic beta", " --linear")
	nom <- ifelse(binary, "logistic", "linear")

	if(file.exists(paste0(out, "/out-", pheno_name, ".assoc.", nom, ".gz")))
	{
		return(NULL)
	}

	cmd <- paste0("plink",
		" --bfile ", geno, 
		" --pheno ", pheno,
		" --pheno-name ", pheno_name,
		flag,
		" --remove ", exclude,
		" --out ", out, "/out-", pheno_name,
		" --covar ", covs,
		" --allow-no-sex",
		" --hide-covar"
	)
	system(cmd)
	unlink(paste0(out, "/out-", pheno_name, ".log"))
	unlink(paste0(out, "/out-", pheno_name, ".nosex"))
	organise_result(paste0(out, "/out-", pheno_name, ".assoc.", nom), fr)
}

organise_result <- function(resfile, fr)
{
	res <- fread(resfile)
	res <- res[,c("CHR", "BP", "TEST"):=NULL]
	res <- left_join(res, fr, by=c("SNP"="ID"))
	index <- res$A1 == res$REF
	res$A2 <- res$ALT
	res$A1[!index] <- res$ALT[!index]
	res$A2[!index] <- res$REF[!index]
	res$ALT_FREQS[index] <- 1 - res$ALT_FREQS[index]
	res$SE <- res$BETA / res$STAT

	res <- select(res, 
		SNP=SNP,
		effect_allele=A1,
		other_allele=A2,
		eaf=ALT_FREQS,
		beta=BETA,
		se=SE,
		pval=P,
		samplesize=NMISS
	)
	write.table(res, file=resfile, row=F, col=T, qu=F)
	system(paste0("gzip ", resfile))
}

##


args <- commandArgs(T)
jid <- as.numeric(args[1])
nphen <- as.numeric(args[2])


##

g <- expand.grid(chunk=1:40, pre=c("binary", "ord", "cont"))
g$filename <- paste0("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/data-", g$pre, "-", g$chunk, "-40.txt")
g$out <- paste0("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/data-", g$pre, "-", g$chunk, "-40-mod.txt")

fn <- list.files("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/mod/")
g <- data.frame(out=paste0("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/mod/", fn), stringsAsFactors=FALSE)
g$binary <- grepl("binary", g$out)

l <- list()
for(i in 1:nrow(g))
{
	message(i)
	if(file.exists(g$out[i]))
	{	
		l[[i]] <- data.frame(
			filename=g$out[i], 
			binary=g$binary[i],
			phen=scan(g$out[i], what=character(), nlines=1)[-c(1:2)]
		)
	}
}

vars <- bind_rows(l)


geno <- "../genotypes/instruments"
fr <- fread(paste0(geno, ".afreq"))
names(fr)[1] <- "chr"
fr <- fr[,c("chr", "OBS_CT"):=NULL]

covs <- "../other_files/covs.txt"
exclude <- "../other_files/exclude.txt"
out <- "../results"

first <- (jid - 1) * nphen + 1
last <- min(nrow(vars), jid * nphen)

vars <- vars[first:last, ]

for(i in 1:nrow(vars))
{
	message(i)
	plink_run(geno, covs, vars$filename[i], vars$phen[i], out, exclude, fr, vars$binary[i])
}


