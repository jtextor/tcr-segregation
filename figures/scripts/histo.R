library( readxl )
library( beeswarm )

source( "tools.R" )


pdf("plots/histo.pdf", width=6.5, height=1.8, pointsize=10, useDingbats=FALSE)

par( mar=c(2,4,2,.2), mgp=normal.plot.mgp, mfrow=c(1,3), bty="l", bty="n", font.main=1, cex=1  )

d <- read_excel( "data/SRBC Histo Ki67 Quantifizierung Rohdaten.xlsx" )

dki67 <- sapply( d[3:5, 2:4], as.numeric )

dki67 <- as.data.frame( dki67 )


nice.beeswarm( dki67, labels=nd34,
	yticks=seq(0,3000,by=1000),
	ylab=expression(paste("Ki67"^"+"," cells" , " / mm"^2)),
       	strictly.positive=TRUE ) 


dgc <- sapply( d[27:29, 2:4], as.numeric )
dgc <- as.data.frame( dgc )

print( dgc )

nice.beeswarm( dgc, yticks=c(.01,.025,.05,.075,.1,.25,.5,.75,1), 
	      ytics.show.each=5,
	      labels=nd34,
	      ylab=expression(paste("GCs / mm"^2)),
	log=TRUE )

dgsize <- sapply( d[32:34, 2:4], as.numeric )
dgsize <- as.data.frame( dgsize ) / 1000

print( dgsize )

nice.beeswarm( dgsize, yticks=c(.1,.25,.5,.75,1,2.5,5,7.5,10,25,50), 
		yticklabels=c("0.1",NA,NA,NA,"1",NA,NA,NA,"10",NA,"50"),
#	      	yticks.show.each=4,
	       labels=nd34,
	      ylab=expression(paste("GC size (1000 mm"^2,")")),
	     log=TRUE )

dev.off()


