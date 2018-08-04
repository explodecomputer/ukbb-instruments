library(dplyr)
library(data.table)

g <- expand.grid(chunk=1:40, pre=c("binary", "ord", "cont"))
g$filename <- paste0("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/data-", g$pre, "-", g$chunk, "-40.txt")
g$out <- paste0("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/data-", g$pre, "-", g$chunk, "-40-mod.txt")

fn <- list.files("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/mod/")
g <- data.frame(out=paste0("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/mod/", fn), stringsAsFactors=FALSE)
g$binary <- grepl("binary", g$out)

l <- list()
for(i in 1:nrow(g))
{
	message(i)
	if(file.exists(g$out[i]))
	{	
		phen <- fread(g$out[i],he=T)
		phen <- phen[,-c(1:2)]
		ncase <- rep(0, ncol(phen))
		ncontrol <- rep(0, ncol(phen))
		for(j in 1:ncol(phen))
		{

			if(g$binary[i])
			{
				ncontrol[j] <- sum(phen[,..j] == 1, na.rm=TRUE)
				ncase[j] <- sum(phen[,..j] == 2, na.rm=TRUE)
			} else {
				ncase[j] <- NA
				ncontrol[j] <- sum(!is.na(phen[,..j]))
			}
		}
		l[[i]] <- data.frame(
			filename=g$out[i], 
			binary=g$binary[i],
			phen=names(phen),
			#phen=scan(g$out[i], what=character(), nlines=1)[-c(1:2)],				
			ncontrol=ncontrol,
			ncase=ncase
		)
	}
}

vars <- bind_rows(l)
save(vars, file="../other_files/phesant_vars.rdata")

