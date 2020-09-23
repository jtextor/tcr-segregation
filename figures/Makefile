.SECONDARY:
.DELETE_ON_ERROR:

all : figures.pdf

figures.pdf : figures.tex plots/basic.pdf plots/clonefreq.pdf plots/vdist.pdf \
	plots/segregation.pdf plots/usage-shared.pdf
	latexmk -use-make -pdf $<

plots/%.pdf : plots-%.R tmp/basic.image tmp/extra-nucs.image tmp/dmat-spleen-naive.tex
	Rscript $<

tmp/dmat-spleen-naive.tex : visualize-segregation.R data/dmat-spleen-naive.txt
	Rscript $<

tmp/extra-nucs.image : plots-vjnuc-prepare.R tmp/overlaps.image
	Rscript $<

tmp/basic.image : plots-basic-prepare.R
	Rscript $<

tmp/overlaps.image : plots-overlaps-prepare.R
	Rscript $<

clean:
	rm -f figures.pdf tmp/* plots/*
