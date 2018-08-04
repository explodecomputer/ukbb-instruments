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

plink_run_quick <- function(geno, covs, pheno, pheno_name, out, exclude, fr, binary)
{
        if(file.exists(paste0(out, "/out-", pheno_name, ".qassoc.gz")))
        {
                return(NULL)
        }
	# Residualise the phenotype
	phen <- fread(pheno, header=TRUE)
	covar <- fread(covs, header=FALSE)
	covar <- subset(covar, select=-c(V1))
	excl <- scan(exclude, what=character())
	nom <- names(covar)[-c(1)] %>% paste(collapse=" + ")
	p <- subset(phen, select=c("FID", pheno_name))
	names(p) <- c("V2", "y")
	p$y <- p$y - 1
	p <- merge(p, covar, by="V2")
	p <- subset(p, !V2 %in% excl)
	fo <- as.formula(paste0("y ~ ", nom))
	if(binary)
	{
		p$res <- residuals(glm(fo, p, family="binomial", na.action=na.exclude))
	} else {
		p$res <- residuals(lm(fo, p, na.action=na.exclude))
	}
	p$res <- scale(p$res)
	phenfile <- paste0(out, "/out-", pheno_name, ".phen")
	write.table(data.frame(p$V2, p$V2, p$res), file=phenfile, row=F, col=F, qu=F)
	cmd <- paste0("plink",
		" --bfile ", geno,
		" --pheno ", phenfile,
		" --assoc",
		" --remove ", exclude,
		" --out ", out, "/out-", pheno_name,
		" --allow-no-sex"
	)
	system(cmd)
    unlink(paste0(out, "/out-", pheno_name, ".log"))
    unlink(paste0(out, "/out-", pheno_name, ".nosex"))
    unlink(paste0(out, "/out-", pheno_name, ".phen"))
    organise_result_quick(paste0(out, "/out-", pheno_name, ".qassoc"), fr)
}

organise_result_quick <- function(resfile, fr)
{
        res <- fread(resfile)
        res <- res[,c("CHR", "BP"):=NULL]
        res <- left_join(res, fr, by=c("SNP"="ID"))
        # index <- res$A1 == res$REF
        # res$A2 <- res$ALT
        res <- select(res,
                SNP=SNP,
                effect_allele=ALT,
                other_allele=REF,
                eaf=ALT_FREQS,
                beta=BETA,
                se=SE,
                pval=P,
                r2=R2,
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

load("../other_files/phesant_vars.rdata")
geno <- "../genotypes/mrbase_instruments_20180612"
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
	# plink_run(geno, covs, vars$filename[i], vars$phen[i], out, exclude, fr, vars$binary[i])
	plink_run_quick(geno, covs, vars$filename[i], vars$phen[i], out, exclude, fr, vars$binary[i])
}


