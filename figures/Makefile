%.tif: %.png
	convert $*.png $*.tif

%.tif: %.pdf
	convert -alpha off -depth 8 -density 300 -compress lzw $*.pdf $*.tif

## should improve Makefile so that we don't rebuild PNG 
figures: figures.R
	#echo "rmarkdown::render('figures.rmd')" | R --slave
	R CMD BATCH figures.R
	./mvfigs

clean:
	rm *.aux *.log
