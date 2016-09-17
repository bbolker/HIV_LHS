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

ltab.full <- cbind(ltab, ltab.het)

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

tvec <- c(1:1000)

I_matFull <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_matFull <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vecFull <- rep(0, n.trial)
peak_matFull <- matrix(NA, nrow = n.trial, ncol = 2)
val_vecFull <- rep(0, n.trial)

if (file.exists("ev_LHS_resFull.rda")) {
	load("ev_LHS_resFull.rda")
} else {
	for (i in 1:100) {
		cat(i)
		HIVpars <- as.HIVvirparlist(ltab.full[i,])
		pp <- expand(HIVpars)
		yini <- calc_yini_full(pp)
		little_r <- get_rval(full_model_rcpp, yini, pp, tvec = c(1:30))
		
		if(little_r < 0.04){
			interval = c(1,10)
		}else{
			interval = c(0.01,1)
		}
		val_vecFull[i] = find_scale(0.04, full_model_rcpp, yini, HIVpars, interval=interval,
															adjpar = "scale_all", full = TRUE, tvec = c(1:300))
		HIVpars_adj <- transform(HIVpars, scale_all = val_vecFull[i])
		pp_adj <- expand(HIVpars_adj)
		r <- rk(unlist(yini), func = full_model_rcpp, parms = pp_adj, times = tvec)
		
		I_matFull[,i] = r[,(ncol(r) - 1)]
		vir_matFull[,i] = r[,ncol(r)]
		eq_vecFull[i] = vir_matFull[length(tvec),i]
		
		peak_matFull[i,1] = which.max(vir_matFull[,i])
		peak_matFull[i,2] = max(vir_matFull[,i])
		
	}
	save("I_matFull", "vir_matFull", "eq_vecFull","peak_matFull", "val_vecFull",file="ev_LHS_resFull.rda")
}
