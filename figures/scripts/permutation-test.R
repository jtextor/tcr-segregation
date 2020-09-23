local({
	load("tmp/fig2.image")
	Mjaccard.naive <<- Mjaccard.naive
	Mjaccard.d3 <<- Mjaccard.d3
	Mjaccard.d4 <<- Mjaccard.d4
})

x <- rep( NA, 6 )

x[1] <- 1

n <- 1

calc.dist <- function(p){
	d <- 0
	for( i in unique(p) ){
		ii <- which(p==i)
		d <- d + M[ii[1],ii[2]]
	}
	R <<- rbind( R, c(p,d) )
}

fill.it <- function( x, i, f ){
	ina <- which(is.na(x))
	if( sum(ina) == 0 ){
		f(x)
	} else {
		x[ina[1]] <- i

		for( j in ina[2:length(ina)] ){
			x[j] <- i
			fill.it( x, i+1, f )
			x[j] <- NA
		}
	
		x[ina[2]] <- NA
	}
}

make.plot <- function(mfile, title,
		    xlab="Assignment of sample to mice",
		    ylab="sum of within-mouse\nsegregated clones"){
	R <<- c()
	if( is.character( mfile ) ){
		M <<- as.matrix(read.table(mfile))
	} else {
		M <<- mfile
	}
	fill.it( rep(NA,ncol(M)), 1, calc.dist )
	cc <- c(2,rep(1,nrow(R)-1))
	cc <- cc[order(R[,7])]
	R <- R[order(R[,7]),]
	barplot( R[,7], main=title, 
		xlab=xlab,
		ylab=ylab,
		col=cc,
		names.arg=1:nrow(R),
		border=NA,
		ylim=c(0,1.3*max(R[,7])) )
}

#par(mfrow=c(2,1))

pdf("plots/permutation-test.pdf",width=5.5,height=5)

par(mfrow=c(2,3), mgp=c(2,.7,0), oma=c(0,1,0,0), font.main=1, cex.main=1)

make.plot("data/dmat-spleen-naive.txt","naive",xlab="")
make.plot("data/dmat-spleen-d3.txt","d3 p.i.",ylab="")
make.plot("data/dmat-spleen-d4.txt","d4 p.i.",xlab="",ylab="")


make.plot(Mjaccard.naive, "naive", ylab="sum within-mouse \n Jaccard indices", xlab="" )
make.plot(Mjaccard.d3, "d3 p.i.", ylab="" )
make.plot(Mjaccard.d4, "d4 p.i.", ylab="", xlab="" )


#make.plot("dmat-srbc.txt","SRBC")

legend("topleft", c("correct", "incorrect"),
	bty="n", fill=c(2,1) )

dev.off()
