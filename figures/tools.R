# Description: parses RTCR output (e.g. "results.tsv") for the R tcR package.
#
# tcR website: http://imminfo.github.io/tcr/

parse.rtcr <- function(.filename)
{
    f <- gzfile(.filename)
    records <- strsplit(tail(readLines(f), n = -1), split = "\t", fixed = T)
    close(f)
    df <- data.frame(Umi.count = NA, Umi.proportion = NA,
                     Read.count = as.integer(sapply(records, function(x)x[1])),
                     CDR3.nucleotide.sequence = sapply(records,
                                                       function(x)x[5]),
                     CDR3.amino.acid.sequence = sapply(records,
                                                       function(x)x[2]),
                     V.gene = sapply(records,
                                     function(x)strsplit(x[3],"*",T)[[1]][1]),
                     D.gene = "",
                     J.gene = sapply(records,
                                     function(x)strsplit(x[4],"*",T)[[1]][1]),
                     V.end = sapply(records, function(x)x[6]),
                     J.start = sapply(records, function(x)x[7]),
                     D5.end = -1,
                     D3.end = -1,
                     VD.insertions = -1,
                     DJ.insertions = -1,
                     Total.insertions = -1,
                     stringsAsFactors = F)
    df$Read.proportion <- df$Read.count / sum(df$Read.count)
    df <- df[, c("Umi.count", "Umi.proportion", "Read.count",
                 "Read.proportion", "CDR3.nucleotide.sequence",
                 "CDR3.amino.acid.sequence", "V.gene", "J.gene", "D.gene",
                 "V.end", "J.start", "D5.end", "D3.end", "VD.insertions",
                 "DJ.insertions", "Total.insertions")]
    cls <- c("as.integer", "as.numeric", "as.integer", "as.numeric",
             "as.character", "as.character", "as.character", "as.character",
             "as.character", "as.integer", "as.integer", "as.integer",
             "as.integer", "as.integer", "as.integer", "as.integer")
    for (i in 1:ncol(df)){
        df[[i]] <- do.call(cls[i], list(df[[i]]))
    }
    df
}

read.rtcr <- function(s="S1"){
	parse.rtcr(paste0("../data/",s,"/results.tsv.gz"))
}


vgenes.order <- c( "V1", "V2" ,  "V3" , "V4" ,   "V5"  ,
     "V12-1", "V12-2", "V13-1", "V13-2", "V13-3", "V14",   "V15",   "V16",  
	"V17",   "V19",     "V20",   "V23",   "V24",   "V26",   "V29", 
	"V30",   "V31"  )


bxax <- function( l ) axis( 1, labels=l, at=seq_along(l), tick=FALSE, line=-1.6, padj=1, mgp=group.plot.mgp )
nd34.axis <- function() bxax(c("naive","d3 p.i.","d4 p.i."))

t.ci <- function( x, conf.level=0.95 ){
	xt <- t.test( x+rnorm(length(x))*0.001, conf.level=conf.level )
	c( xt$conf.int[1], xt$estimate, xt$conf.int[2] )
}

nice.beeswarm <- function( M, labels=colnames(M), yticks, yticklabels=as.character(yticks),
			  yticks.show.each=NULL,
			  ylab, conf.level=0.95, pch=19,
			 strictly.positive=FALSE, log=FALSE, ... ){
	if( max(sapply(M,max)) > max(yticks) || min(sapply(M,min)) < min(yticks) ){
		stop("Values would be truncated!")
	}
	if( log ){
		M <- log(M) 
		ytickst <- log( yticks )
	} else {
		ytickst <- yticks
	}
	Mc <- as.data.frame(lapply( M, t.ci ))
	if( strictly.positive == TRUE ){
		Mc[1,Mc[1,]<0] <- 0
	}
	bxplot( Mc, probs=c(0,0.5,1), ylim=range(ytickst), width=.5, at=seq(1,ncol(Mc)),
	      labels="", xaxt="n", yaxt="n", ylab=ylab,  ... )
	beeswarm( M, labels="", pch=pch, add=TRUE, ... )
	#	 yaxt="n", xaxt="n",
#		pch=19,
#		ylab=ylab, ...)


	bxax( labels )
	if( !is.null( yticks.show.each ) ){
		yticklabels[ -seq(1,length(yticklabels),by=yticks.show.each) ] <- NA
	}	
	axis(2, at=ytickst, labels=yticklabels)
}

nd34 <-  c("naive","d3 p.i.","d4 p.i.")

compare.distributions <- function( xl, col=rep("gray33",length(xl)),
		xticks=seq(0,15,by=5), xticks.each=3, yticks.each=5,
		ylim=c(0,20) ){
	mb <- function(x, x0, w=4){
		x <- sum(as.double(names(x))*x)/100
		arrows( x0, x, x0-.2*diff(range(xticks)), x, length=.02, col="red" ) 
	}

	col <- rep( col, length(xl) / length(col) )

	x <- xl[[1]]
	
	grps <- as.character( seq( ylim[1], ylim[2] ) )

	bc <- barplot(x[grps],xlim=range(xticks), xaxt="n", yaxt="n", 
		col=col[1],
		horiz=TRUE, border=NA, ylab="",
		xlab="", main="shared")
	lbls <- as.character( xticks )
	lbls[-seq(1,length(lbls),by=xticks.each)] <- NA

	axis( 1, at=xticks, labels=lbls)

	mb(x,max(xticks))

	grpls <- grps
	grpls[-seq(1,length(grpls),by=yticks.each)] <- NA

	axis( 2, at=bc, labels=grpls, tick=FALSE )

	i<-2
	for( x in xl[-1] ){
		barplot(x[grps],horiz=TRUE,xlim=range(xticks), col=col[i],
			border=NA,main=c("unique","shared")[i%%2+1], xaxt="n", yaxt="n")
		axis( 1, at=xticks, labels=lbls)
		mb(x,max(xticks))
		i<-i+1
	}
}

mouse.colors <- c("lightgray","darkgray","black") #c("red","blue","gold")
mouse.colors.tikz <- c("black!33","black!66","black")

betweenmice.pch <- 1

normal.plot.mgp <- c(2,.7,0)
group.plot.mgp <- c(2,1,0)


