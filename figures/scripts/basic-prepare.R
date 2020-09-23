
library( tcR )

merge.nuc.clones <- function(x) tapply( x$Read.count, list(x$"CDR3.amino.acid.sequence"), FUN=sum )

olap <- function( i1, i2 ){
	x <- tapply( naive.data[[i1]]$Read.count,
		list(naive.data[[i1]]$"CDR3.amino.acid.sequence"), FUN=sum )
	y <- tapply( naive.data[[i2]]$Read.count,
		list(naive.data[[i2]]$"CDR3.amino.acid.sequence"), FUN=sum )
	on <- intersect( names(x), names(y) )
	length(on) / (length(x) + length(y))
}

diver <- function(x) chao1(tapply(x$Read.count, list(x$"CDR3.amino.acid.sequence"), FUN=sum))


source("tools.R")
naive.data <- lapply(paste0("S",1:6),read.rtcr)
repseq.stats(naive.data)
rs <- cloneset.stats(naive.data)
d <- sapply( naive.data, diver )

save.image( file="tmp/basic.image" )

