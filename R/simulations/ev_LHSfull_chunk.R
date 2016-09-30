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

t0 <- proc.time()
argvals <- commandArgs(trailingOnly=TRUE) ## chunk number: full run is (20)x50
## indexed from 0
batch_num <- as.numeric(argvals[1])
batch_size <- as.numeric(argvals[2])
cat(batch_num,batch_size,"\n")

start_sim <- batch_num*batch_size+1  
end_sim <- (batch_num+1)*batch_size
fn <- paste0("ev_LHSfull_",batch_num,"_v2.rda")
cat(start_sim,end_sim,fn,"\n")

n.trial <- end_sim-start_sim+1 ## overwrite source()d value

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

tvec <- 1:2000

I_matFull <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_matFull <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vecFull <- rep(0, n.trial)
peak_matFull <- matrix(NA, nrow = n.trial, ncol = 2)
val_vecFull <- rep(0, n.trial)

mrow <- 1
for (i in start_sim:end_sim) {
    cat(i,proc.time()-t0,"\n")
    r <- try(
        { HIVpars <- as.HIVvirparlist(ltab[i,])
            pp <- expand(HIVpars)
            yini <- calc_yini_full(pp)
            little_r <- get_rval(full_model_rcpp, yini, pp)
            
            if (little_r < 0.04) {
                interval = c(1,10)
            } else {
                interval = c(0.01,1)
            }
            val_vecFull[mrow] = find_scale(0.04, full_model_rcpp, yini,
                                        HIVpars, interval=interval,
                                        adjpar = "scale_all", full = TRUE)
            HIVpars_adj <- transform(HIVpars, scale_all = val_vecFull[mrow])
            pp_adj <- expand(HIVpars_adj)
            rk(unlist(yini), func = full_model_rcpp, parms = pp_adj,
               times = tvec, hmax = 0.5)
            })
    if (!inherits(r,"try-error")) {
        I_matFull[,mrow] = r[,(ncol(r) - 1)]
        vir_matFull[,mrow] = r[,ncol(r)]
        eq_vecFull[mrow] = vir_matFull[length(tvec),mrow]
    
        peak_matFull[mrow,1] = which.max(vir_matFull[,mrow])
        peak_matFull[mrow,2] = max(vir_matFull[,mrow])

        save("I_matFull", "vir_matFull", "eq_vecFull","peak_matFull", "val_vecFull",file=fn)
    }
    mrow <- mrow+1
}
save("I_matFull", "vir_matFull", "eq_vecFull","peak_matFull", "val_vecFull",file=fn)

