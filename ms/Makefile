%.pdf: %.rmd
	echo "library(rmarkdown); render(\"$*.rmd\")" | R --slave

%.pdf: %.Rnw
	echo "library(knitr); knit2pdf(\"$*.Rnw\")" | R --slave

%.pdf: %.tex
	texi2dvi -p $*.tex

%.tex: %.Rnw
	echo "library(knitr); knit(\"$*.Rnw\")" | R --slave

ms%.pdf: ms$*.Rnw ms_parms.tex
	echo "library(knitr); knit2pdf(\"ms$*.Rnw\")" | R --slave


compare: plos_latex.tex
	../latexdiff/latexdiff plos_latex_old.txt plos_latex.tex >plos_latex_cmp.tex
	texi2dvi --batch plos_latex_cmp.tex
clean:
	rm -f ms.html ms.tex .#* #*
