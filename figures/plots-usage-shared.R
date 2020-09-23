load("tmp/fig2.image")

source("tools.R")

library( beeswarm )

pdf("plots/fig2.pdf", width=2.3, height=2, pointsize=10, useDingbats=FALSE )

par( mar=c(2,4,1.5,0.5), mgp=normal.plot.mgp, bty="n"  )

jacc.plot <- function( M ){
	d.same <- 100*c( M[1,2], M[3,4], M[5,6] )
	d.diff <- 100*c( M[1,3:6], M[2,3:6], M[3,5:6], M[4,5:6] )

	nice.beeswarm( list( d.same, d.diff ), yticks=seq(0,12,by=2), 
		labels=c("within\nmice","between\nmice"), ylab="Jaccard index (%)", pch=19,
		pwcol=c(mouse.colors,rep("black",12)),
		pwpch=c(rep(19,3),rep(betweenmice.pch,12)) )
}

jacc.plot( Mjaccard.naive )
jacc.plot( Mjaccard.d3 )
jacc.plot( Mjaccard.d4 )

dev.off()

prop.ci <- function( x, conf.level=0.95 ){
	x <- log( x / ( 1-x ) )
	xt <- t.test( x+rnorm(length(x))*0.001, conf.level=conf.level )
	r <- c( xt$conf.int[1], xt$estimate, xt$conf.int[2] )
	exp(r) / (exp(r)+1)
}

colPropCi <- function( M ){
	100*apply( M/100, 2, prop.ci )
}

usage.plots <- function( xr.v, xr.j, .legend=NULL, group=FALSE, each=3,
		       col=c("gray33","gray66"), ylab="usage (% of reads)", level="within") {

	par( fig = c(0,.65,0,1) )

	xr.v <- xr.v[,vgenes.order]

	if( group ){
		f <- rep(1:2,each=each)
		xr.v.err <- rbind( colPropCi(xr.v[f==1,]), colPropCi(xr.v[f==2,]) )
		xr.v <- rbind( colMeans(xr.v[f==1,]), colMeans(xr.v[f==2,]) )


		xr.j.err <- rbind( colPropCi(xr.j[f==1,]), colPropCi(xr.j[f==2,]) )
		xr.j <- rbind( colMeans(xr.j[f==1,]), colMeans(xr.j[f==2,]) )

		each <- 1 
	}

	bc <- barplot( xr.v, beside=T, las=2, ylab=ylab, border=NA, col=rep(col,each=each),
		ylim=c(0,25) )
	if( group ){
		segments( bc, matrix(xr.v.err[c(1,4),],nrow=1), bc, matrix(xr.v.err[c(3,6),],nrow=1) )
	}
	if( !is.null(.legend) ){
		legend( "topright", .legend, col=col, pch=16, bty="n" ) 
	}
	if( !is.null( level ) ){
		mtext( paste(level," mice"), 3 )
	}

	par( fig = c(.65,1,0,1), new=T )

	barplot( xr.j, beside=T, las=2, ylab=ylab, border=NA, col=rep(col, each=each),
		ylim=c(0,25) ) 
	

	if( group ){
		segments( bc, matrix(xr.j.err[c(1,4),],nrow=1), bc, matrix(xr.j.err[c(3,6),],nrow=1) )
	}

	if( !is.null( level ) ){
		mtext( paste(level," mice"), 3 )
	}
}


pdf("plots/usage-shared.pdf",width=6.5,height=2,pointsize=10, useDingbats=FALSE)

par( mar=c(4,4,1,0.5), mgp=c(mgp=c(2,.7,0)), oma=c(0,0,0,0) )

usage.plots(xr.v.wi,xr.j.wi,.legend=c("shared","unique"), group=TRUE)
#mtext( "within mice", 3, outer=TRUE )

usage.plots(xr.v,xr.j, level="between", group=TRUE)
#mtext( "between mice", 3, outer=TRUE )


usage.plots(xr.v.nd3, xr.j.nd3, .legend=c("naive","d3 p.i."), group=TRUE, each=6, level=NULL )

dev.off()

pdf("plots/reads-shared.pdf",width=2,height=2,pointsize=10, useDingbats=FALSE)

par( mar=c(2,4,1.5,0.5), mgp=c(mgp=c(2,.7,0)), bty="n" )

#plot( density( log( tapply( nsw$"Number of reads", list(nsw$"Amino acid sequence"), FUN=sum ) ) ) )
#lines( density( log( tapply( nuw$"Number of reads", list(nuw$"Amino acid sequence"), FUN=sum ) ) ) )

boxplot( readcounts.wi,
	names=c("shared","unique"), 
	ylab=expression(paste(log[10]," number of reads") ), outline=FALSE, xaxt="n", ylim=c(0,4) )
axis( 1, lty=0, labels=c("shared","unique"), at=c(1,2) )
title( "within mice", font.main=1, line=0.5, cex.main=1 )

boxplot( readcounts.bet,
	names=c("shared","unique"), xaxt="n",
	ylab=expression(paste(log[10]," number of reads") ), outline=FALSE, ylim=c(0,4) )
axis( 1, lty=0, labels=c("shared","unique"), at=c(1,2) )



title( "between mice", font.main=1, line=0.5, cex.main=1 )


dev.off()


