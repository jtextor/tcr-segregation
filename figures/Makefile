.SECONDARY:
.DELETE_ON_ERROR:

all : figures.pdf

figures.pdf : figures.tex plots/basic.pdf plots/clonefreq.pdf plots/vdist.pdf \
	plots/segregation.pdf plots/usage-shared.pdf plots/permutation-test.pdf \
	plots/fig4.pdf plots/dge.pdf plots/histo.pdf
	latexmk -use-make -pdf $<

plots/%.pdf : scripts/%.R tmp/basic.image tmp/extra-nucs.image tmp/dmat-spleen-naive.tex \
	tmp/dge.image
	Rscript $<

tmp/dmat-spleen-naive.tex : scripts/visualize-segregation.R data/dmat-spleen-naive.txt tmp/overlaps.image
	Rscript $<

tmp/extra-nucs.image : scripts/vjnuc-prepare.R tmp/overlaps.image
	Rscript $<

tmp/dge.image : scripts/dge-prepare.R tmp/overlaps.image
	Rscript $<

tmp/%.image : scripts/%-prepare.R
	Rscript $<

clean:
	rm -f figures.pdf tmp/* plots/*
