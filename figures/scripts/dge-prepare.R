
library(data.table)
	
fr <- function(f,seq.col="Amino acid sequence", read.col="Number of reads"){ 
	g <- fread(f,data.table=FALSE)[,c(seq.col, read.col)]
	aggregate( g[,2], list(aaseq=g[,1]), sum )
}
	

if( TRUE ){
	  
	rd0 <- Reduce( function(...) merge(...,all=TRUE,by="aaseq"), 
		lapply( paste0("../data/B",1:3,"/results.tsv.gz"), fr ) )
	
	rd3 <- Reduce( function(...) merge(...,all=TRUE,by="aaseq"), 
		lapply( paste0("../data/B3",1:3,"/results.tsv.gz"), fr ) )
	
	rd4 <- Reduce( function(...) merge(...,all=TRUE,by="aaseq"), 
		lapply( paste0("../data/B4",1:3,"/results.tsv.gz"), fr ) )
	
	colnames(rd0) <- c("aaseq","N_1","N_2","N_3")
	colnames(rd3) <- c("aaseq","D3_1","D3_2","D3_3")
	colnames(rd4) <- c("aaseq","D4_1","D4_2","D4_3")
	
	#r.all <- merge( rd0, merge( rd3, rd4, all=TRUE ), all=TRUE )
	
	rr <- merge( rd0, rd3, all=TRUE )
	chck <- rowSums(apply(rr,2,is.na))
	rr <- rr[chck < 5,]
	
	rrr <- merge(rr,rd4, all=TRUE)
	chck <- rowSums(apply(rrr,2,is.na))
	rrr <- rrr[chck < 8,]
	
	rrr[is.na(rrr)] <- 0
	write.csv( rrr, file="tmp/merged-blood.csv" )
}

if( TRUE ){

	cat("read d0 data ..." )
	rd0 <- Reduce( function(...) merge(...,all=TRUE,by="aaseq"), 
		lapply( paste0("../data/S",1:6,"/results.tsv.gz"), fr ) )
	colnames(rd0) <- c("aaseq","N_1","N_2","N_3","N_4","N_5","N_6")
	cat("\n")

	cat("read d3 data ..." )
	rd3 <- Reduce( function(...) merge(...,all=TRUE,by="aaseq"), 
		lapply( paste0("../data/S3",1:6,"/results.tsv.gz"), fr ) )
	colnames(rd3) <- c("aaseq","D3_1","D3_2","D3_3","D3_4","D3_5","D3_6")
	cat("\n")
	
	cat("read d4 data ..." )
	rd4 <- Reduce( function(...) merge(...,all=TRUE,by="aaseq"), 
		lapply( paste0("../data/S4",1:6,"/results.tsv.gz"), fr ) )
	colnames(rd4) <- c("aaseq","D4_1","D4_2","D4_3","D4_4","D4_5","D4_6")
	cat("\n")

	rr <- merge( rd0, rd3, all=TRUE )
	rr <- merge( rr, rd4, all=TRUE )

	rr[is.na(rr)] <- 0

	rr <- rr[rowSums(rr[,-1]>0)>1,]
	write.csv( rr, file="tmp/merged-spleen.csv" )
}

dgel <- function(f="blood"){
	library( edgeR )
	r <- as.matrix(read.csv(paste0("tmp/merged-",f,".csv"),row.names="aaseq")[,-1])
	r <- r[rowSums(r>0)>=ncol(r)/2,]

	#r <- r[,c(1:3,7:9)]
	
	condition <- factor(
		c( rep("naive",ncol(r)/3), rep("ag", 2*ncol(r)/3 ) ) ) #, rep("d4", 3) ) )
	
	y <- DGEList( counts=r, group=condition )
	y <- calcNormFactors(y)
	design <- model.matrix( ~condition )
	y <- estimateDisp( y, design )

	fit <- glmFit( y, design, prior.count=1 )

	tst <- glmLRT( fit, coef=2 )

	list(y, tst)
}

lB <- dgel()
lS <- dgel("spleen")

save(lB, lS, file="tmp/dge.image")

