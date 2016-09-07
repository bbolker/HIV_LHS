library("deSolve")
library("reshape2")
library("ggplot2"); theme_set(theme_bw())
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

I_mat16 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_mat16 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vec16 <- rep(0, n.trial)
peak_mat16 <- matrix(NA, nrow = n.trial, ncol = 2)
val_vec16 <- rep(0, n.trial)

if (file.exists("ev_LHS_res16.rda")) {
  load("ev_LHS_res16.rda")
} else {
  for (i in 1:n.trial) {
    cat(i)
    HIVpars <- as.HIVvirparlist(ltab[i,])
    yini <- calc_yini(HIVpars)
    little_r <- get_rval(gfun2(HIVpars),yini, HIVpars)

    if(little_r < 0.04){
      
      interval = c(1,10)
    }else{
      interval = c(0.01,1)
    }
    val_vec16[i] = find_scale(0.04, gfun2, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vec16[i])
    r <- lsoda(unlist(yini), func = gfun2(HIVpars_adj), parms = HIVpars_adj, times = tvec)

    I_mat16[,i] = r[,(ncol(r) - 1)]      
    vir_mat16[,i] = r[,ncol(r)]
    eq_vec16[i] = vir_mat16[length(tvec),i]
    
    peak_mat16[i,1] = which.max(vir_mat16[,i])
    peak_mat16[i,2] = max(vir_mat16[,i])
    
  }
  save("I_mat16","vir_mat16", "eq_vec16","peak_mat16", "val_vec16",file="ev_LHS_res16.rda")
}
