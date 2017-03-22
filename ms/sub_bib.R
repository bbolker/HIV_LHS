txt <- readLines("plos_latex.tex")
bibtxt <- readLines("plos_latex.bbl")
target <- grep("\\bibliography{virulence}",txt,fixed=TRUE)
if (length(target)==0) stop("can't find bib line")
out_txt <- c(txt[1:(target-1)],
             bibtxt,
             txt[target:length(txt)])
writeLines("plos_latex_out.tex")
