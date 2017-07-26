library(data.table)

check_binary <- function(x)
{
	temp <- head(table(x))
	stopifnot(length(temp) == 2)
	print(temp)
	if(0 %in% names(temp) & 1 %in% names(temp))
	{
		message("adding")
		x <- x + 1
	} else if(1 %in% names(temp) & 2 %in% names(temp)) {
		message("fine")
	} else {
		message("subbing")
		x[x==names(temp)[1]] <- 9999999
		x[x==names(temp)[2]] <- 9999998
		x[x==9999999] <- 0
		x[x==9999998] <- 1
	}
	return(x)
}

setwd("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/")

# f <- list.files(pattern="data*")
# f <- f[!grepl("unord", f)]

g <- expand.grid(chunk=1:40, pre=c("binary", "ord", "cont"))
g$filename <- paste0("data-", g$pre, "-", g$chunk, "-40.txt")
g$out <- paste0("data-", g$pre, "-", g$chunk, "-40-mod.txt")

mapping <- fread("../../data.7445.csv", he=T, sep=",")

for(i in 1:nrow(g))
{
	message(i)
	a <- fread(g$filename[i], header=TRUE, sep=",")
	if(ncol(a) > 1)
	{
		index <- match(a$userID, mapping$app15825)
		stopifnot(all(a$userID == mapping$app15825[index]))
		a$userID <- mapping$app8786[index]
		for(j in 2:ncol(a))
		{
			if(g$pre[i] == "binary")
			{
				a[,j] <- check_binary(a[,..j])

			}
		}
		temp <- data.frame(FID=a$userID, IID=a$userID)
		a <- cbind(temp, a[,c("userID"):=NULL])
		write.table(a, file=g$out[i], row=FALSE, col=TRUE, qu=FALSE)
	}
}


