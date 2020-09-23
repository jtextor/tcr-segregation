library( beeswarm )
source( "tools.R" )

p <- function(f="dmat-spleen-naive.txt"){
	Dm <- as.matrix(read.table(paste0("data/",f)))
	seg.1.2 <- c(Dm[c(1,2),c(3,4)])
	seg.1.3 <- c(Dm[c(1,2),c(5,6)])
	seg.2.3 <- c(Dm[c(3,4),c(5,6)])
	c(seg.1.2,seg.1.3,seg.2.3)
}

p.wi <- function(f="dmat-spleen-naive.txt"){
	Dm <- as.matrix(read.table(paste0("data/",f)))
	Dm[cbind(c(1,3,5),c(2,4,6))]
}

pdf("plots/segregation.pdf", width=2, height=2, pointsize=10, useDingbats=FALSE)

par( mar=c(2,4,1,.5), mgp=normal.plot.mgp, bty="l", bty="n"  )

nice.beeswarm( list( p.wi(), p() ), labels=c("within\nmice","between\nmice"), ylab="segregated clones",
	pwcol=c(mouse.colors,rep("black",12)), pwpch=c(rep(19,3),rep(betweenmice.pch,12)),
	yticks=seq(0,40,by=10) )

nice.beeswarm( list( p.wi(), p.wi("dmat-spleen-d3.txt") ), 
	labels=nd34[1:2],
	pwcol=rep(mouse.colors,2), pwpch=c(rep(19,3),rep(17,3)),ylab="segregated clones",
	yticks=seq(0,300,by=100))
mtext( "within mice", 3, line=0 )

nice.beeswarm( list( p(), p("dmat-spleen-d3.txt") ),
	labels=nd34[1:2], 
	pch=betweenmice.pch,
	ylab="segregated clones", yticks=seq(0,300,by=100) )
mtext( "between mice", 3, line=0 )

dev.off()


pdf("plots/segregation2.pdf", width=3.25, height=1.5, pointsize=10, useDingbats=FALSE )

par( mar=c(2,4,1,.5), mgp=normal.plot.mgp, bty="l", bty="n" )

nice.beeswarm( list( p.wi(), p.wi("dmat-spleen-d3.txt"), p.wi("dmat-spleen-d4.txt") ), 
	labels=nd34,
	pwcol=rep(mouse.colors,3), pwpch=c(rep(19,3),rep(17,3),rep(15,3)),
	ylab="segregated clones", yticks=seq(0,300,by=100)  )
mtext( "within mice", 3, line=0 )


nice.beeswarm( list( p(), p("dmat-spleen-d3.txt"), p("dmat-spleen-d4.txt") ), 
	pch=betweenmice.pch, 
	labels=nd34,
	ylab="segregated clones", yticks=seq(0,300,by=100)  )
mtext( "between mice", 3, line=0 )


dev.off()


