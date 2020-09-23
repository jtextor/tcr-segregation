

library( data.table )

naive.data <- lapply(paste0("S",1:6), function(x) fread( paste0("../data/",x,"/results.tsv.gz"),data.table=FALSE ))
d3.data <- lapply(paste0("S3",1:6), function(x) fread( paste0("../data/",x,"/results.tsv.gz"),data.table=FALSE ))
d4.data <- lapply(paste0("S4",1:6), function(x) fread( paste0("../data/",x,"/results.tsv.gz"),data.table=FALSE ))

olap <- function( i1, i2 ){
	x <- tapply( current.tcr.data[[i1]]$"Number of reads", list(current.tcr.data[[i1]]$"Amino acid sequence"), FUN=sum )
	y <- tapply( current.tcr.data[[i2]]$"Number of reads", list(current.tcr.data[[i2]]$"Amino acid sequence"), FUN=sum )
	on <- intersect( names(x), names(y) )
	#sum( c(x,y)[on] ) / sum(c(x,y))
	length(on) / (length(x) + length(y))
}



current.tcr.data <- naive.data
M <- matrix( 0, 6, 6 )
apply( combn( 1:6, 2 ), 2, function(x) { M[x[1],x[2]] <<- M[x[2],x[1]] <<- olap( x[1], x[2] ) } )
Mjaccard.naive <- M

current.tcr.data <- d3.data
M <- matrix( 0, 6, 6 )
apply( combn( 1:6, 2 ), 2, function(x) { M[x[1],x[2]] <<- M[x[2],x[1]] <<- olap( x[1], x[2] ) } )
Mjaccard.d3 <- M


current.tcr.data <- d4.data
M <- matrix( 0, 6, 6 )
apply( combn( 1:6, 2 ), 2, function(x) { M[x[1],x[2]] <<- M[x[2],x[1]] <<- olap( x[1], x[2] ) } )
Mjaccard.d4 <- M



shared <- function( a, b ){
	da <- do.call( rbind, naive.data[a] )
	db <- do.call( rbind, naive.data[b] )
	xn <- intersect( da$"Amino acid sequence",
		db$"Amino acid sequence" )
	rbind( da[ da$"Amino acid sequence" %in% xn, ], 
		db[ db$"Amino acid sequence" %in% xn, ] )
}


unique <- function( a, b ){
	da <- do.call( rbind, naive.data[a] )
	db <- do.call( rbind, naive.data[b] )
	xn <- intersect( da$"Amino acid sequence",
		db$"Amino acid sequence" )
	rbind( da[ ! da$"Amino acid sequence" %in% xn, ], 
		db[ ! db$"Amino acid sequence" %in% xn, ] )
}


v.usage <- function(x) tapply( x$"Number of reads", list(x$"V gene"), FUN=function(y) 100*sum(y)/sum(x$"Number of reads") )
j.usage <- function(x) tapply( x$"Number of reads", list(x$"J gene"), FUN=function(y) 100*sum(y)/sum(x$"Number of reads") )

usage.forplot <- function(  .shared, .unique, .usagefun ) {
	xr <- do.call( rbind, lapply( c(.shared, .unique), .usagefun ) )
	colnames(xr) <- substring( colnames(xr),4)
	colnames(xr) <- gsub( "\\*.*", "", colnames(xr) )
	xr
}

xr.v.nd3 <- usage.forplot( naive.data, d3.data, v.usage )
xr.j.nd3 <- usage.forplot( naive.data, d3.data, j.usage )

naive.shared.wi <- lapply( 1:3, function(x) shared( 2*x-1, 2*x ) ) 
naive.unique.wi <- lapply( 1:3, function(x) unique( 2*x-1, 2*x ) ) 

xr.v.wi <- usage.forplot( naive.shared.wi , naive.unique.wi, v.usage )
xr.j.wi <- usage.forplot( naive.shared.wi , naive.unique.wi, j.usage )


naive.shared <- list( shared( 1:2, 3:4 ), shared( 1:2, 5:6 ), shared( 3:4, 5:6 ) )
naive.unique <- list( unique( 1:2, 3:4 ), unique( 1:2, 5:6 ), unique( 3:4, 5:6 ) )

xr.v <- usage.forplot( naive.shared , naive.unique, v.usage )
xr.j <- usage.forplot( naive.shared , naive.unique, j.usage )



nsw <- do.call( rbind, naive.shared.wi ) 
nuw <- do.call( rbind, naive.unique.wi ) 

readcounts.wi <- list( x=unname( log10( tapply( nsw$"Number of reads", list(nsw$"Amino acid sequence"), FUN=sum ) ) ),
	y=unname( ( log10( tapply( nuw$"Number of reads", list(nuw$"Amino acid sequence"), FUN=sum ) ) ) ) )
rm( nsw )
rm( nuw )

nsb <- do.call( rbind, naive.shared ) 
nub <- do.call( rbind, naive.unique ) 

readcounts.bet <- list( x=unname( log10( tapply( nsb$"Number of reads", list(nsb$"Amino acid sequence"), FUN=sum ) ) ),
	y=unname( ( log10( tapply( nub$"Number of reads", list(nub$"Amino acid sequence"), FUN=sum ) ) ) ) )

rm( nsb )
rm( nub )

rm( naive.data )
rm( d3.data )
rm( d4.data )

extract.aaseq <- function(x) x[,c("Number of reads","Junction nucleotide sequence","V gene end position","J gene start position")]

naive.unique.aaseq <- lapply( naive.unique, extract.aaseq )
naive.unique.wi.aaseq <- lapply( naive.unique.wi, extract.aaseq )

naive.shared.aaseq <- lapply( naive.shared, extract.aaseq )
naive.shared.wi.aaseq <- lapply( naive.shared.wi, extract.aaseq )

save( naive.unique.aaseq, naive.shared.aaseq,
     naive.unique.wi.aaseq, naive.shared.wi.aaseq, file="tmp/overlaps.image" )

rm( naive.unique.aaseq, naive.shared.aaseq,
     naive.unique.wi.aaseq, naive.shared.wi.aaseq )

rm( naive.shared.wi )
rm( naive.unique.wi )

rm( naive.shared )
rm( naive.unique )

rm( current.tcr.data )


save.image( file="tmp/fig2.image" )


