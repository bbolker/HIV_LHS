%.tif: %.png
	convert $*.png $*.tif

%.tif: %.pdf
	convert -alpha off -depth 8 -density 600 -compress lzw $*.pdf $*.tif

## rasterizes, but does pretty well?
%.eps: %.pdf
	inkscape --export-area-page --export-eps=$*.eps $*.pdf


tikz-figs: partnerSIR.tex diagrams.tex diagram_geom.tex
	. ./mktikzfigs

tikz-comb: comb_diagrams.tex partnerSIR.tex diagrams.tex diagram_geom.tex
	pdflatex comb_diagrams
	pdfcrop comb_diagrams
	pdftops -eps comb_diagrams-crop.pdf
	mv comb_diagrams-crop.eps fig2.eps
	mv comb_diagrams-crop.pdf fig2.pdf

## should improve Makefile so that we don't rebuild PNG 
figures: figures.R
	R CMD BATCH --vanilla figures.R
##	./mvfigs

clean:
	rm *.aux *.log
