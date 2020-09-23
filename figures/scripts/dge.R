load("tmp/dge.image")

source("tools.R") 
library( edgeR )
library( beeswarm )

plot.new()

pdf("plots/dge.pdf", width=6.5, height=3, pointsize=10, useDingbats=FALSE )

par( mfcol=c(2,3), mar=c(2,4,2,0.2), cex.main=1.0, 
    mgp=c(2,.7,0), font.main=1, bty="n", cex=1 )

recs <- rownames( topTags( lB[[2]], 3 ) )


for( rec in recs ){
	nice.beeswarm( as.data.frame(matrix(cpm(lB[[1]],log=T,prior.count=1)[rec,],ncol=3)),
		yticks=c(0,5,10,15,20),
		strictly.positive=TRUE,
		pwcol=rep(rep(mouse.colors,each=1),3),
		pwpch=c(rep(19,3),rep(17,3),rep(15,3)),	
		labels=nd34,
		main=paste(rec,"\n(blood)"),
		ylab="" 
	)


	if( rec == recs[1] ){
		title( ylab=expression( paste( log[2],' counts / 10'^6 ) ) ) 
	}

	nice.beeswarm(as.data.frame(matrix(cpm(lS[[1]],log=T,prior.count=1)[rec,],ncol=3)),
		yticks=seq(0,15,by=5),
		strictly.positive=TRUE,
		pwcol=rep(rep(mouse.colors,each=2),3),
		pwpch=c(rep(19,6),rep(17,6),rep(15,6)),	
		labels=nd34,
		main=paste(rec,"\n(spleen)"),
		ylab=""
	)

	if( rec == recs[1] ){
		title( ylab=expression( paste( log[2],' counts / 10'^6 ) ) ) 
	}
}

dev.off()

