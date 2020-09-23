load("tmp/basic.image")

source("tools.R")

pdf("plots/clonefreq.pdf", width=1.7, height=3, pointsize=10, useDingbats=FALSE )

par(bty="n", mar=c(4,3,1,0), mgp=normal.plot.mgp)
 #,mfrow=c(1,3))

clonefreq <- function(x) {
	x <- tapply( x$Read.count, list(x$"CDR3.amino.acid.sequence"), FUN=sum )
	t1 <- table( x ) 
	list( x = as.integer(names(t1)), y = t1 )
}

plot( clonefreq( rbind(naive.data[[1]],naive.data[[2]]) ), 
	pch=19, cex=.3, col=mouse.colors[1],log="xy", 
     xlab="reads in clone", ylab="", ylim=c(1,10000) )
points( clonefreq( rbind(naive.data[[3]],naive.data[[4]]) ), pch=19, cex=.3, col=mouse.colors[2] ) 
points( clonefreq( rbind(naive.data[[5]],naive.data[[6]]) ), pch=19, cex=.3, col=mouse.colors[3] ) 
legend( "topright", c("1","2", "3"), title="mouse", col=mouse.colors, pch=19, bty="n" )

mtext("number of clones",2,2.2)

dev.off()

