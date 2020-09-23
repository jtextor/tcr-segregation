load("tmp/basic.image")

source("tools.R")

library( beeswarm )

pdf("plots/basic.pdf", width=1.8, height=2.4, pointsize=10, useDingbats=FALSE )

par(bty="n", mar=c(2,4,1,0), mgp=c(2,.7,0), yaxs='r' ) #,mfrow=c(1,3))

xl <- function(){
	axis( 1, at=1:3, tick=FALSE, line=-1 )
	title( xlab="mouse", line=1 )
}

rs[,"Sum.reads"] <- rs[,"Sum.reads"]/1000

beeswarm(rs[,"Sum.reads"]~rep(1:3,each=2), xlab="", 
	xaxt="n",
	ylab="",
	ylim=c(0,max(2000,max(rs[,"Sum.reads"])))#,pwpch=rep(c(1,19),3)
	, pch=19, col=mouse.colors)
xl()

mtext( "mapped reads", 2, line=3 )
mtext( "(x1000)", 2, line=1.8 )

rs[,"#Aminoacid clonotypes"] <- rs[,"#Aminoacid clonotypes"] / 1000

beeswarm(rs[,"#Aminoacid clonotypes"]~rep(1:3,each=2), xlab="", 
	#ylab=expression(paste("unique CDR3",beta," AA sequences")),
	xaxt="n",
	ylab="",
	ylim=c(0,max(80,rs[,"#Aminoacid clonotypes"])),
	#,pwpch=rep(c(1,19),3)
	pch=19, col=mouse.colors)
xl()

mtext( expression(paste("unique CDR3",beta,"")), 2, line=2.7 )
mtext( " AA sequences (x1000)", 2, line=1.8 )

d[1,] <- d[1,]/1000

beeswarm(d[1,]~rep(1:3,each=2), 
	xlab="", xaxt="n",
	ylab="",
	ylim=c(0,80),
	#pwpch=rep(c(1,19),3),
	pch=19, col=mouse.colors)

xl()

mtext( expression(paste("corrected unique CDR3",beta,"")), 2, line=2.7 )
mtext( " AA sequences (x1000)", 2, line=1.8 )


dev.off()


