source("tools.R")

mk <- function( f = "dmat-spleen-naive", fout=paste0("tmp/",f,".tex") ){

	if( is.character( f ) ){
		m <- as.matrix( read.table( paste0("data/",f,".txt") ) )
	} else {
		m <- f
	}

	m <- 1+log(m+1)

#print(m)
#print( diag(m) )
#print(m[row(m) != col(m)])

	m[row(m) != col(m)] <- scale(m[row(m) != col(m)])

#print( min(m) )

	m <- m - min(m) +1 # +min(m[row(m) != col(m)])

	#print(m)

	if( !is.null( fout ) ){
		sink( fout )
	}

	for( i in 1:nrow(m) ){
		for( j in 1:ncol(m) ){
			if( i > j ){
				col <- ""
				if( j %% 2 == 1 && i == j+1 ){
					if( j == 1 ){
						col <- paste0(',densely dashed,draw=',mouse.colors.tikz[1])
					} 
					if( j == 3 ){
						col <- paste0(',densely dashed,draw=',mouse.colors.tikz[2])
					}
					if( j == 5 ){
						col <- paste0(',densely dashed,draw=',mouse.colors.tikz[3])
					}
				}
				cat( "\\draw [line width=",(m[i,j]),"pt",col,"] (x", i, ") -- (x", j, ");\n", sep="" )
			}
		}
	}

	if( !is.null( fout ) ){
		sink()
	}
}

mk( "dmat-spleen-naive" )
mk( "dmat-spleen-d3" )
mk( "dmat-spleen-d4" )

load( "tmp/fig2.image" )
mk( M, "tmp/dmat-spleen-jaccard-naive.tex" )

