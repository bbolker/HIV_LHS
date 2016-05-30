hivlhs_eeid2016.tex: hivlhs_eeid2016.Rnw
	echo "knitr::knit('hivlhs_eeid2016.Rnw')" | R --slave
	sed -i.bak -e 's/\\documentclass{beamer}/\\documentclass[final,hyperref={pdfpagelabels=false}]{beamer}/' hivlhs_eeid2016.tex

hivlhs_eeid2016.pdf: hivlhs_eeid2016.tex
	texi2dvi -p hivlhs_eeid2016.tex

figures: figure/univ.R figure/tradeoff.R figure/
