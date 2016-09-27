library("Rcpp")
library("deSolve")
library("reshape2")
library("ggplot2"); theme_set(theme_bw())
sourceCpp("../fullModel.cpp")
source("../fullModel.R")
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")

batch_num <- BATCHNUM  ## chunk number: full run is (20)x50
## indexed from 0
start_sim <- batch_num*20+1  
end_sim <- (batch_num+1)*20
fn <- paste0("ev_LHSfull_",batch_num,".rda")
cat(start_sim,end_sim,fn,"\n")

n.trial <- end_sim-start_sim+1 ## overwrite source()d value

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

tvec <- 1:2000

I_matFull <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_matFull <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vecFull <- rep(0, n.trial)
peak_matFull <- matrix(NA, nrow = n.trial, ncol = 2)
val_vecFull <- rep(0, n.trial)


for (i in start_sim:end_sim) {
    cat(i)
    HIVpars <- as.HIVvirparlist(ltab[i,])
    pp <- expand(HIVpars)
    yini <- calc_yini_full(pp)
    little_r <- get_rval(full_model_rcpp, yini, pp)
    
    if(little_r < 0.04){
        interval = c(1,10)
    }else{
        interval = c(0.01,1)
    }
    val_vecFull[i] = find_scale(0.04, full_model_rcpp, yini, HIVpars, interval=interval,
                                adjpar = "scale_all", full = TRUE)
    HIVpars_adj <- transform(HIVpars, scale_all = val_vecFull[i])
    pp_adj <- expand(HIVpars_adj)
    r <- rk(unlist(yini), func = full_model_rcpp, parms = pp_adj, times = tvec)
    
    I_matFull[,i] = r[,(ncol(r) - 1)]
    vir_matFull[,i] = r[,ncol(r)]
    eq_vecFull[i] = vir_matFull[length(tvec),i]
    
    peak_matFull[i,1] = which.max(vir_matFull[,i])
    peak_matFull[i,2] = max(vir_matFull[,i])

    save("I_matFull", "vir_matFull", "eq_vecFull","peak_matFull", "val_vecFull",file=fn)

}
save("I_matFull", "vir_matFull", "eq_vecFull","peak_matFull", "val_vecFull",file=fn)

