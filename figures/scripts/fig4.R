library( data.table )

p <- function( a=1, b=2, pad=.2 ){

foundinB.A <- names( head(naive.data[[a]],20 ) ) %in% names( naive.data[[b]] )
foundinA.B <- names( head(naive.data[[b]],20 ) ) %in% names( naive.data[[a]] )

plot( pad+head(naive.data[[a]],20 ), ylim=c(-1.25,1.225), pch=1+18*foundinB.A, bty="n", 
	xlab="", ylab="",
	xaxt="n", yaxt="n" )
points( -pad-head(naive.data[[b]],20 ), pch=1+18*foundinA.B )

common.top20 <- intersect( names( head(naive.data[[a]],20 ) ), names( head(naive.data[[b]],20 ) ) )
common.A <- match( common.top20, names(naive.data[[a]]) )
common.B <- match( common.top20, names(naive.data[[b]]) )

segments( common.A, pad+naive.data[[a]][common.A], common.B, -pad-naive.data[[b]][common.B] )

abline( h=pad, col="gray" )
abline( h=-pad, col="gray" )


}

ax <- function( pad=.2, ymax=1, by=1 ){
	xt <- c( -seq(0,ymax,by=by)-pad, seq(0,ymax,by=by) + pad )
	#axis( 2, at=xt, labels=abs(rep(seq(0,ymax,by=by),2)) )
	axis( 2, at=seq(0,ymax,by=by)+pad, labels=paste0(seq(0,ymax,by=by),"%") )
	title( ylab="TCZ 1", adj=0.8, line=2 )
	axis( 2, at=-seq(0,ymax,by=by)-pad, labels=paste0(seq(0,ymax,by=by),"%") )
	title( ylab="TCZ 2", adj=0.2, line=2 )
}

plt <- function(){
p(1,2)
ax()
legend( "bottomright", c("present in other TCZ","absent in other TCZ"), pch=c(19,1), bty="n")
title( "mouse 1", font.main=1, cex.main=1, line=-1 )
p(3,4)
ax()
title( "mouse 2", font.main=1, cex.main=1, line=-1 )
title( ylab="share of sample", line=3.5 )
p(5,6)
title( "mouse 3", font.main=1, cex.main=1, line=-1 )
axis( 1, at=c(1,5,10,15,20) )
ax()
title( xlab="clone index", line=2, outer=TRUE )
}


getdata <- function(pre="S"){
	naive.data <<- lapply(paste0(pre,1:6), function(x) fread( paste0("../data/",x,"/results.tsv.gz") ))

	naive.data <<- lapply( naive.data,
		function( x) tapply( x$"Number of reads", list(x$"Amino acid sequence"), FUN=sum ) )

	naive.data <<- lapply( naive.data, function(x) 100 * sort( x, decreasing=T ) / sum(x) )

}

pdf("plots/fig4.pdf", width=3.5,height=5.5, pointsize=10, useDingbats=FALSE )

par( mfrow=c(3,1), mar=c(2,4.5,0,0), oma=c(0,0,2,0),
	mgp=c(mgp=c(2,.7,0)), cex=1 )


getdata("S")
plt()
par( oma=c(0,4.5,2,0) )
title( main="naive animal", font.main=1, cex.main=1, outer=TRUE )
par( oma=c(0,0,2,0) )

getdata("S3")
plt()
par( oma=c(0,4.5,2,0) )
title( main="day 3 post infection", font.main=1, cex.main=1, outer=TRUE )
par( oma=c(0,0,2,0) )


getdata("S4")
plt()
par( oma=c(0,4.5,2,0) )
title( main="day 4 post infection", font.main=1, cex.main=1, outer=TRUE )
par( oma=c(0,0,2,0) )


dev.off()

