library(data.table)

check_binary <- function(x)
{
	temp <- head(table(x))
	stopifnot(length(temp) == 2)
	# print(temp)
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


split_by_col <- function(x, fn, nc=3)
{
	n <- ncol(x)-2
	if(n <= nc)
	{
		x$userID <- paste(x$userID, x$userID)
		names(x)[1] <- "FID IID"
		write.table(x, file=fn, row=FALSE, col=TRUE, qu=FALSE)
	} else {
		x <- as.data.frame(x)
		temp <- data.frame(FID=x$userID, IID=x$userID)
		x <- subset(x, select=-c(userID))
		l <- list()
		nsplit <- ceiling(n/nc)
		for(i in 1:nsplit)
		{
			message(i)
			i1 <- (i-1) * nc + 1
			i2 <- min(i * nc, n)
			fn1 <- gsub(".txt", paste0(".", i, ".txt"), fn)
			x1 <- cbind(temp, x[,i1:i2])
			write.table(x1, file=fn1, row=FALSE, col=TRUE, qu=FALSE)
		}
	}
}

setwd("~/data/bb/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant/")

# f <- list.files(pattern="data*")
# f <- f[!grepl("unord", f)]

g <- expand.grid(chunk=1:40, pre=c("catord", "cont", "binary"))
g$filename <- paste0("data-", g$pre, "-", g$chunk, "-40.txt")
g$out <- paste0("mod/data-", g$pre, "-", g$chunk, "-40.txt")

mapping <- fread("../../data.7445.csv", he=T, sep=",")

for(i in 112:nrow(g))
{
	message(i, " of ", nrow(g))
	a <- fread(g$filename[i], header=TRUE, sep=",")
	dim(a)
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

		split_by_col(a, g$out[i], nc=100)
	}
}


