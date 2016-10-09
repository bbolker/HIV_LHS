## go to directory with individual batchfile outputs ...
## setwd("../../simdata/ev_LHSfull")

tag <- "_v3"  ## batch version identifier
nrun <- function(x) {
   r <- readLines(x)
   sum(grepl("([0-9.]+ ){5}",r))
}
Routfiles <- list.files(pattern=paste0(tag,"\\.Rout"))
if (length(Routfiles)>0) {
    outruns <- sapply(Routfiles,nrun)
    sum(outruns)
}

nsucc <- function(x) {
   load(x)
   sum(!apply(is.na(I_matFull),2,all))
}

rdafiles <- list.files(pattern=paste0(tag,"\\.rda"))
if (length(rdafiles)==0) {
    warning("no .rda files; are you in the right working directory?")
} else {
    table(ss <- sapply(rdafiles,nsucc))

    sum(ss)  ## 213

    get_all <- function(x) {
        ## cat(x,"\n")
        L <- load(x)
        r <- lapply(L,function(x) get(x,envir=parent.frame(2)))
        names(r) <- L
        return(r)
    }

    comb_fun <- function(x,nruns=20) {
        vals <- lapply(all_list,"[[",x)
        vals <- vals[!sapply(vals,function(x) all(is.na(x)))]
        if (!is.null(dd <- dim(all_list[[1]][[x]]))) {
            if (dd[2]==nruns) {
                res <- do.call(cbind,vals)
            } else {
                res <- do.call(rbind,vals)
            }
        } else res <- do.call(c,vals)
        return(res)
    }

    all_list <- lapply(rdafiles,get_all)
    c_names <- names(all_list[[1]])
    all_comb <- lapply(c_names,comb_fun)
    names(all_comb) <- c_names
    save("all_comb",file=paste0("../ev_LHSfull_comb",tag,".rda"))
}

if (FALSE) {
    ## setwd("..") ## back up to main simdata dir
    load(paste0("ev_LHSfull_comb",tag,".rda"))
    ## quick look at results ...
    pcol <- adjustcolor("black",alpha=0.1)
    matplot(all_comb$I_matFull,type="l",lty=1,col=pcol, log = "y")
    
    ## make sure initial growth rate is consistent even when equilibrium prevalence is close to subset limits
    matplot(all_comb$I_matFull[,all_comb$I_matFull[2000,] < 2e-2], type = "l", lty = 1, col = pcol, log = "y")
    ## These seem OK although they seem a bit off. Should we get rid of them? or keep them? Re-adjust little r just for these?
    
    matplot(all_comb$vir_matFull,type="l",lty=1,col=pcol,
            ylim=c(2.5,7))
    
    sum(apply(all_comb$vir_matFull,
          2,function(x)
        !all(is.na(x)) & min(x)>2.5  & max(x)<7))
    ## 785 (out of 792)
    plot(density(all_comb$eq_vecFull))
    plot(density(na.omit(all_comb$peak_matFull[,1])))
    plot(density(na.omit(all_comb$peak_matFull[,2])))
    plot(density(all_comb$val_vecFull))
}
