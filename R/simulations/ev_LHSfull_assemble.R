tag <- "_v2"  ## batch version identifier
nrun <- function(x) {
   r <- readLines(x)
   sum(grepl("([0-9.]+ ){5}",r))
}
Routfiles <- list.files(pattern=paste0(tag,"\\.Rout"))
outruns <- sapply(Routfiles,nrun)
sum(outruns)

nsucc <-
function(x) {
   load(x)
   sum(!apply(is.na(I_matFull),2,all))
}

rdafiles <- list.files(pattern=paste0(tag,"\\.rda"))
table(ss <- sapply(rdafiles,nsucc))

sum(ss)  ## 143

get_all <- function(x) {
  cat(x,"\n")
  L <- load(x)
  return(setNames(lapply(L,get),L))
}

all_list <- lapply(rdafiles,get_all)
lapply(names(all_list[[1]]),
    function(x) {
       if (!is.null(dd <- dim(all_list[[1]][[x]]))) {
          do.call(rbind,lapply(all_list,"[[",x))
       }
     })
