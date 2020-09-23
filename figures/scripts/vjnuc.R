
source("tools.R")
load("tmp/extra-nucs.image" )

pdf("plots/vjnuc.pdf", width=3.3, height=2.5, pointsize=10, useDingbats=FALSE)

par( mfrow=c(1,4), font.main=1,cex.main=1, mar=c(4,1,1,0), mgp=c(2,.7,0), oma=c(0,4,1,1), cex=1)

compare.distributions( vj.distances, col=c("gray33","gray66") )

par( oma=c(4,4,1,1) )
title( ylab="nucleotides between \n V and J segments", outer=TRUE, line=1 )
title( xlab="% of reads", outer=TRUE )
mtext( "within mice", adj=0.2, outer=TRUE )
mtext( "between mice", adj=0.9, outer=TRUE )



par( oma=c(0,4,1,1) )
compare.distributions( extra.nucs, col=c("gray33","gray66"), 
	xticks=seq(0,30,by=10), xticks.each=3, yticks.each=2,
	ylim=c(0,10) )
par( oma=c(4,4,1,1) )
title( ylab="non-D nucleotides between \n V and J segments", outer=TRUE, line=1 )
title( xlab="% of reads", outer=TRUE )
mtext( "within mice", adj=0.2, outer=TRUE )
mtext( "between mice", adj=0.9, outer=TRUE )

dev.off()
