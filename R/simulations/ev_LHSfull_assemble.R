## go to directory with individual batchfile outputs ...
## setwd("../../simdata/ev_LHSfull")

tag <- "_v2"  ## batch version identifier
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
    save("all_comb",file="../ev_LHSfull_comb.rda")
}

if (FALSE) {
    ## setwd("..") ## back up to main simdata dir
    load("ev_LHSfull_comb.rda")
    ## quick look at results ...
    pcol <- adjustcolor("black",alpha=0.1)
    matplot(all_comb$I_matFull,type="l",lty=1,col=pcol)
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
    
    ##little_r
    
    min(all_comb$I_matFull[2000,], na.rm = TRUE)
    
    Imat <- all_comb$I_matFull[,!is.na(all_comb$I_matFull[1,])]
    
    little_r <- unlist(apply(Imat, 2, function(x){
    	df <- data.frame(t = 1:2000, I = x)
    	l <- try(lm(log(I)~t, data = df, subset = I > 1e-3 & I < 1e-2))
    	ifelse(!inherits(l,"try-error"), coef(l)[[2]], NA) 
    }))
    
    plot(density(little_r, na.rm = TRUE))
    matplot(Imat[,which(little_r < 0.05)], type = "l", lty = 1, col = pcol, log = "y")
    length(which(little_r < 0.05))
    length(which(little_r > 0.05))
}
