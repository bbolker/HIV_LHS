library("deSolve")
library("reshape2")
library("ggplot2"); theme_set(theme_classic())
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")
library("gridExtra")

get_rval <- function(g, yini, pp, plot.it=FALSE,
                     tvec = c(1:500),
                     lims=c(0.01,0.05),
                     verbose=FALSE) {
  
    start <- unlist(yini)
    
    r <- rk(y=start,
                 times=tvec,
                 func=g,
                 parms=pp, method = "ode45")
        
    Itot <- r[,(ncol(r)-1)]
    
    if(max(Itot) < 5e-2){
      return(0)
    }else{
      
      dd <- data.frame(tvec=tvec,Itot=Itot)
      
      dsub <- subset(dd, Itot>lims[1] & Itot<lims[2])
      mm <- try(lm(log(Itot)~tvec,data=dsub))
      if (plot.it) {
        plot(log(Itot)~tvec,data=dsub)
        abline(mm,col=2)
      }
      cc <- coef(mm)[2]
      return(cc)
    }
}

HIVpars.shirreff <- transform(HIVpars.shirreff, ini_V = 3)

tvec = seq(from =1, to =600, by = 0.1)

I_matS <- matrix(NA, nrow = length(tvec), ncol = 3)
vir_matS <- matrix(NA, nrow = length(tvec), ncol = 3)
val_vecS <- rep(NA, 3)

ini_r <- c(0.021, 0.042, 0.084)

if (file.exists("ev_LHS_resS.rda")) {
  load("ev_LHS_resS.rda")
} else {
  for(i in 1:3){
    cat(i)
    HIVpars <- transform(HIVpars.shirreff, Ini_I = 0.001)
    yini <- calc_yini3(HIVpars)
    
    little_r <- get_rval(gfun5(HIVpars),yini, HIVpars)
    
    if(little_r < ini_r[i]){
      interval = c(1,2)
    }else{
      interval = c(0.01,1)
    }
    
    val_vecS[i] = find_scale(ini_r[i], gfun5, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vecS[i])
    r <- rk(unlist(yini), func = gfun5(HIVpars), parms = HIVpars_adj, times = tvec, method = "ode45")
    
    I_matS[,i] = r[,(ncol(r) - 1)]
    vir_matS[,i] = r[,ncol(r)]
  }
  save("I_matS", "vir_matS", "val_vecS",file="ev_LHS_resS.rda")
}

I_matS2 <- matrix(NA, nrow = length(tvec), ncol = 3)
vir_matS2 <- matrix(NA, nrow = length(tvec), ncol = 3)
val_vecS2 <- rep(NA, 3)

iniI_vec <- c(0.001, 0.0001, 0.00001)

if (file.exists("ev_LHS_resS2.rda")) {
  load("ev_LHS_resS2.rda")
} else {
  for(i in 1:3){
    cat(i)
    HIVpars <- transform(HIVpars.shirreff, Ini_I = iniI_vec[i])
    yini <- calc_yini3(HIVpars)
    
    little_r <- get_rval(gfun5(HIVpars),yini, HIVpars)
    
    if(little_r < 0.042){
      interval = c(1,2)
    }else{
      interval = c(0.01,1)
    }
    
    val_vecS2[i] = find_scale(0.042, gfun5, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vecS2[i])
    r <- lsoda(unlist(yini), func = gfun5(HIVpars), parms = HIVpars_adj, times = tvec)
    
    I_matS2[,i] = r[,(ncol(r) - 1)]
    vir_matS2[,i] = r[,ncol(r)]
  }
  save("I_matS2", "vir_matS2", "val_vecS2",file="ev_LHS_resS2.rda")
}

I_matS3 <- matrix(NA, nrow = length(tvec), ncol = 3)
vir_matS3 <- matrix(NA, nrow = length(tvec), ncol = 3)
val_vecS3 <- rep(NA, 3)

iniI_vir <- c(2.5, 3, 3.5)

if (file.exists("ev_LHS_resS3.rda")) {
  load("ev_LHS_resS3.rda")
} else {
  for(i in 1:3){
    cat(i)
    HIVpars <- transform(HIVpars.shirreff, Ini_I = 0.001, ini_V = iniI_vir[i])
    yini <- calc_yini3(HIVpars)
    
    little_r <- get_rval(gfun5(HIVpars),yini, HIVpars)
    
    if(little_r < 0.042){
      interval = c(1,2)
    }else{
      interval = c(0.01,1)
    }
    
    val_vecS3[i] = find_scale(0.042, gfun5, yini, HIVpars, interval=interval, adjpar = "scale_all")
    HIVpars_adj <- transform(HIVpars, scale_all = val_vecS3[i])
    r <- lsoda(unlist(yini), func = gfun5(HIVpars), parms = HIVpars_adj, times = tvec)
    
    I_matS3[,i] = r[,(ncol(r) - 1)]
    vir_matS3[,i] = r[,ncol(r)]
  }
  save("I_matS3", "vir_matS3", "val_vecS3",file="ev_LHS_resS3.rda")
}