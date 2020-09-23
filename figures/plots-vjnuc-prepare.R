library( Biostrings )
load( "tmp/overlaps.image" )

d1 <- DNAString("GGGACAGGGGGC")
d2 <- DNAString("GGGACTGGGGGGGC")

M <- nucleotideSubstitutionMatrix( match=5, mismatch=-9, baseOnly=TRUE )



xs <- function(x, N=1000){
	xs <- substr( x$"Junction nucleotide sequence", x$"V gene end position"+1, x$"J gene start position"-1 )
	x <- sample( xs, N, replace=TRUE, prob=x$"Number of reads" )
	empty <- x==""
	x <- x[!empty]
	xd <- sapply( x, DNAString )

	a1 <- pairwiseAlignment( pattern = xd, 
	       gapOpening=-Inf, gapExtension=-Inf,
	       substitutionMatrix=M,  
	       subject = d1, type="local"  )

	a2 <- pairwiseAlignment( pattern = xd, 
	       gapOpening=-Inf, gapExtension=-Inf,
	       substitutionMatrix=M,  
	       subject = d2, type="local"  )

	r <- rep(0,N)
	r[!empty] <- unname( apply( rbind( nchar(x)-nmatch(a1), nchar(x)-nmatch(a2) ), 2, min ) )

	r

	#list( x, a2 )
}

extra.nucs <- list( 
	xs( do.call( rbind, naive.shared.wi.aaseq ), 1000 ),		
	xs( do.call( rbind, naive.unique.wi.aaseq ), 1000 ),

	xs( do.call( rbind, naive.shared.aaseq ), 1000 ),
	xs( do.call( rbind, naive.unique.aaseq ), 1000 )
)

extra.nucs <- lapply( extra.nucs, function(x){ x <- table(x); 100*x/sum(x) } )


prep <- function(d){
	d <- do.call( rbind, d )
	d$d <- d$"J gene start position" - d$"V gene end position" - 1
	d <- tapply( d$"Number of reads", d$d, FUN=sum )
	d <- 100*d/sum(d)
}

x1 <- prep( naive.shared.wi.aaseq )
x2 <- prep( naive.unique.wi.aaseq )
x3 <- prep( naive.shared.aaseq )
x4 <- prep( naive.unique.aaseq )

vj.distances <- list( x1, x2, x3, x4 )

save( extra.nucs, vj.distances, file="tmp/extra-nucs.image" )
