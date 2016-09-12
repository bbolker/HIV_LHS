library("deSolve")
library("reshape2")
library("ggplot2"); theme_set(theme_bw())
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

I_mat20 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_mat20 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vec20 <- rep(0, n.trial)
peak_mat20 <- matrix(NA, nrow = n.trial, ncol = 2)
val_vec20 <- rep(0, n.trial)

if (file.exists("ev_LHS_res20.rda")) {
  load("ev_LHS_res20.rda")
} else {
  for (i in 1:n.trial) {
    cat(i)
    HIVpars <- as.HIVvirparlist(ltab[i,])
    yini <- calc_yini3(HIVpars)
    little_r <- get_rval(gfun_random(HIVpars),yini, HIVpars)
    
    if(little_r < 0.04){
      
      interval = c(1,5)
    }else{
      interval = c(0.01,1)
    }
    val_vec20[i] = find_scale(0.04, gfun_random, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vec20[i])
    r <- lsoda(unlist(yini), func = gfun_random(HIVpars_adj), parms = HIVpars_adj, times = tvec)
    
    I_mat20[,i] = r[,(ncol(r) - 1)]
    vir_mat20[,i] = r[,ncol(r)]
    eq_vec20[i] = vir_mat20[length(tvec),i]
    
    peak_mat20[i,1] = which.max(vir_mat20[,i])
    peak_mat20[i,2] = max(vir_mat20[,i])
    
  }
  save("I_mat20", "vir_mat20", "eq_vec20","peak_mat20", "val_vec20", file="ev_LHS_res20.rda")
}
