%.pdf: %.rmd
	echo "library(rmarkdown); render(\"$*.rmd\")" | R --slave

%.pdf: %.Rnw
	echo "library(knitr); knit2pdf(\"$*.Rnw\")" | R --slave

%.tex: %.Rnw
	echo "library(knitr); knit(\"$*.Rnw\")" | R --slave

ms%.pdf: ms$*.Rnw ms_parms.tex
	echo "library(knitr); knit2pdf(\"ms$*.Rnw\")" | R --slave
