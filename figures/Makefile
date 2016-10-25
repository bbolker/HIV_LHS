%.tif: %.png
	convert $*.png $*.tif

%.tif: %.pdf
	convert -alpha off -depth 8 -density 300 -compress lzw $*.pdf $*.tif

## should improve Makefile so that we don't rebuild PNG 
figures: figures.rmd
	echo "rmarkdown::render('figures.rmd')" | R --slave
	./mvfigs
