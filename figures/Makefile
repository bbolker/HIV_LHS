%.tif: %.png
	convert $*.png $*.tif

## should improve Makefile so that we don't rebuild PNG 
figures: figures.rmd
	echo "rmarkdown::render('figures.rmd')" | R --slave
	./mvfigs
