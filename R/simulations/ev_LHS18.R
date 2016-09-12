library("deSolve")
library("reshape2")
library("ggplot2"); theme_set(theme_bw())
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")

HIVpars.skeleton <- transform(HIVpars.skeleton, ini_V = 3)

I_mat18 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
vir_mat18 <- matrix(NA, nrow = length(tvec), ncol = n.trial)
eq_vec18 <- rep(0, n.trial)
peak_mat18 <- matrix(NA, nrow = n.trial, ncol = 2)
val_vec18 <- rep(0, n.trial)

if (file.exists("ev_LHS_res18.rda")) {
  load("ev_LHS_res18.rda")
} else {
  for (i in 1:n.trial) {
    cat(i)
    HIVpars <- as.HIVvirparlist(ltab[i,])
    yini <- calc_yini2(HIVpars)
    little_r <- get_rval(gfun4(HIVpars),yini, HIVpars)

    if(little_r < 0.04){
      
      interval = c(1,10)
    }else{
      interval = c(0.01,1)
    }
    val_vec18[i] = find_scale(0.04, gfun4, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vec18[i])
    r <- lsoda(unlist(yini), func = gfun4(HIVpars_adj), parms = HIVpars_adj, times = tvec)
    
    I_mat18[,i] = r[,(ncol(r) - 1)]
    vir_mat18[,i] = r[,ncol(r)]
    eq_vec18[i] = vir_mat18[length(tvec),i]
    
    peak_mat18[i,1] = which.max(vir_mat18[,i])
    peak_mat18[i,2] = max(vir_mat18[,i])
    
  }
  save("I_mat18", "vir_mat18", "eq_vec18","peak_mat18", "val_vec18",file="ev_LHS_res18.rda")
}
