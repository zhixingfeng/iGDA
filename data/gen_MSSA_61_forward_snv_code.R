library(seqinr)
x <- read.table('MSSA_61_forward.m5', header=FALSE,  as.is=TRUE)

idx.snv <- list()
base.snv <- list()
for (i in 1:nrow(x)){
	cur.shift <- rep(0,nchar(x[i,18]))
	cur.shift[s2c(x[i,19])=='-'] <- 1
	cur.shift <- cumsum(cur.shift)
	
	cur.idx <- which(s2c(x[i,17])!=s2c(x[i,19]) & s2c(x[i,17])!='-' & s2c(x[i,19])!='-')
	idx.snv[[i]] <- as.numeric(x[i,8]) + cur.idx  - cur.shift[cur.idx]
	base.snv[[i]] <- s2c(x[i,17])[cur.idx]
}

names(idx.snv) <- x[,1]

system('rm -f MSSA_61_forward_snv_code.benchmark')
for (i in 1:length(idx.snv))
{
	cat(paste(paste(idx.snv[[i]],base.snv[[i]],sep=','),collapse='\t'), file='MSSA_61_forward_snv_code.benchmark', append=TRUE)
	cat('\t\n', file='MSSA_61_forward_snv_code.benchmark', append=TRUE)
}


