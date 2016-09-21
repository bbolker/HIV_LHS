library("deSolve")
library("reshape2")
library("ggplot2"); theme_set(theme_bw())
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

I_mat17 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_mat17 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vec17 <- rep(0, n.trial)
peak_mat17 <- matrix(NA, nrow = n.trial, ncol = 2)
val_vec17 <- rep(0, n.trial)

if (file.exists("ev_LHS_res17.rda")) {
  load("ev_LHS_res17.rda")
} else {
  for (i in 1:n.trial) {
    cat(i)
    HIVpars <- as.HIVvirparlist(ltab[i,])
    yini <- calc_yini2(HIVpars)
    little_r <- get_rval(gfun3(HIVpars),yini, HIVpars)
    
    if(little_r < 0.04){
      
      interval = c(1,10)
    }else{
      interval = c(0.01,1)
    }
    val_vec17[i] = find_scale(0.04, gfun3, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vec17[i])
    r <- rk(unlist(yini), func = gfun3(HIVpars_adj), parms = HIVpars_adj, times = tvec)
    
    I_mat17[,i] = r[,(ncol(r) - 1)]
    vir_mat17[,i] = r[,ncol(r)]
    eq_vec17[i] = vir_mat17[length(tvec),i]
    
    peak_mat17[i,1] = which.max(vir_mat17[,i])
    peak_mat17[i,2] = max(vir_mat17[,i])
    
  }
  save("I_mat17","vir_mat17", "eq_vec17","peak_mat17", "val_vec17",file="ev_LHS_res17.rda")
}
