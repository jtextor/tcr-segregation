
load("tmp/basic.image")

source("tools.R")

pdf("plots/vdist.pdf", width=3.3, height=3, pointsize=10 )

par( bty="l", mar=c(4,3,1,0), mgp=normal.plot.mgp)

v.usage <- function(x) tapply( x$Read.count, list(x$V.gene), FUN=function(y) 100*sum(y)/sum(x$Read.count) )
j.usage <- function(x) tapply( x$Read.count, list(x$J.gene), FUN=function(y) 100*sum(y)/sum(x$Read.count) )

xr <- do.call( rbind, lapply( naive.data, v.usage ) )
colnames(xr) <- substring( colnames(xr),4)

xr <- xr[,vgenes.order]

barplot( xr, beside=T, las=2, col=rep(mouse.colors,each=2), border=NA,
	ylab="usage (% of reads)", ylim=c(0,25) ) 

dev.off()

pdf("plots/jdist.pdf", width=2, height=3, pointsize=10 )

par( bty="l", mar=c(4,3,1,0), mgp=normal.plot.mgp)

xr <- do.call( rbind, lapply( naive.data, j.usage ) )
colnames(xr) <- substring( colnames(xr),4)
xr <- xr[,sort(colnames(xr))]



barplot( xr, beside=T, las=2, col=rep(mouse.colors,each=2), 
	border=NA, ylab="usage (% of reads)", ylim=c(0,25) )
dev.off()

