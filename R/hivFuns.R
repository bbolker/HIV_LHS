expand <- function(x, ...) {
    UseMethod("expand")
}
expand.HIVvirparlist <- function(x,...) {
    ## mutate() doesn't work this deeply ...
    x <- within(x,
       {
          min.alpha <- alphaDist[["min"]]
          max.alpha <- alphaDist[["max"]]
          d.alpha <- alphaDist[["delta"]]
          alpha <- seq(alphaDist[["min"]],alphaDist[["max"]],by=d.alpha)
          n.alpha <- length(alpha)
          
          rho <- scale_all * scale_c * rho_base
          c_mean <- scale_all * scale_c * c_mean_base
          
          v <- 10^alpha
          ## duration <- (Dmax * D50^h)/(D50^h + v^h)
          Duration2 <- hill(v,Dmax,D50,h)
          lam <- 1/(Duration1 + Duration2 + Duration3)
          lammat <- matrix(rep(lam,n.alpha), n.alpha, n.alpha)
          lammat_dis <- lammat
          diag(lammat_dis) <- 2 * diag(lammat_dis)
          lammat_adj <- lammat + t(lammat)
          ## Beta2 <- Bmax * v^a / (v^a + K^a)
          Beta2 <- hill(K,Bmax,v,a)
          Beta <- scale_all * scale* (Beta1 * Duration1 + Beta2 * Duration2 + Beta3 * Duration3) * lam
          c_e <- Beta * c_e_ratio
          c_u <- Beta * c_u_ratio
          
          ##Risk distribution
          gamma.shape <- 1/kappa
          gamma.scale <- mu/gamma.shape
          bin.risk <-  qgamma(seq(0.01, 0.99, length.out = n.risk+1), shape = gamma.shape, scale = gamma.scale)
          
          p.risk <- rep(1/n.risk, n.risk)
          r.risk <- (bin.risk[-1] + bin.risk[-(n.risk+1)])/2
          
          ##Heterogeneity
          alpha2 <- rep(alpha, n.risk)
          
          cc_mat <- outer(r.risk * c_mean, r.risk * c_mean, "+")/2
          rho2 <- r.risk * rho
          c2 <- r.risk * c_mean
          rho2.I <- rep(rho2, each = n.alpha)
          c2.I <- rep(c2, each = n.alpha)
          cc_mat.I <- outer(c2, c2.I, "+")/2
          cc_mat.II <- outer(c2.I, c2.I, "+")/2
          
          #Infection
          p <- matrix(NA, n.alpha, n.alpha)
          for(i in 1:n.alpha){
          	x <- min.alpha + (i-1)*d.alpha
          	denom <- 
          		diff(pnorm(c(min.alpha-d.alpha/2,max.alpha+d.alpha/2),
          							 mean=x, sd=Vm))
          	for (j in 1:n.alpha){
          		k <- min.alpha + (j-1)*d.alpha
          		p[i,j] = diff(pnorm(k + c(-1,1)*d.alpha/2, mean = x, sd = Vm))/denom
          	}
          }
          
          beta.SI <- matrix(rep(Beta, n.risk*n.risk), nrow = n.risk, byrow = TRUE)
          
          #extracouple
          c_u_ratio2 <- r.risk * c_u_ratio
          c_e_ratio2 <- r.risk * c_e_ratio
          c_u_ratio2.I <- rep(c_u_ratio2, each = n.alpha)
          c_e_ratio2.I <- rep(c_e_ratio2, each = n.alpha)
          c_u2 <- outer(r.risk, c_u)
          r.risk.couple <- outer(r.risk, r.risk, "+")
          
          lam2 <- rep(lam, n.risk)
          lammat_dis2 <- matrix(rep(lam2, n.alpha*n.risk), n.alpha * n.risk, n.alpha * n.risk)
          diag(lammat_dis2) <- 2 * diag(lammat_dis2)
          lammat_adj2 <- outer(lam2, lam2, "+")
          
       })
    return(x)
}

hill <- function(x,a,b,p) {
    a*b^p/(b^p+x^p)
}

##' get r value within specified bounds
##'
##' @param g
##' @param yini initial state
##' @param pp parameters
##' @param plot.it plot linear fit?
##' @param tvec
##' @param lims
##' @param verbose verbose (debugging) output?
get_rval <- function(g, yini, pp, plot.it=FALSE,
                     tvec = c(1:500),
                     lims=c(1e-3,1e-2),
                     verbose=FALSE) {
  
    start <- unlist(yini)
    
    r <- rk(y=start,
                 times=tvec,
                 func=g,
                 parms=pp, hmax = 0.3)
        
    Itot <- r[,(ncol(r)-1)]
    
    if(max(Itot) < 1e-2){
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
      
      if(is.na(cc)){
      	return(3)
      }
      
      return(cc)
    }
}
## get_rval(gfun3(HIVpars1), yini4, HIVpars1) 0.01358363

get_rval2 <- function(val, g, yini, basepar,
                      adjpar="scale",
                      verbose=FALSE,
											full = FALSE, ...) {
    pp <- basepar
    pp[[adjpar]] <- val
    
    if(full){
    	pp2 <- expand(pp)
    	r <- get_rval(g, yini, pp2,...)
    }else{
    	r <- get_rval(g(pp), yini, pp,...)
    }
    if (verbose) cat(adjpar,r,"\n")
    return(r)
}

## get_rval2(5, gfun2, yini3, HIVpars1) 0.022
## get_rval2(1.2, gfun3, yini4, HIVpars1) 0.02492347
## get_rval2(2.1, gfun, yini1, HIVpars1, mutation = TRUE) 0.01887
## get_rval2(4, gfun, yini1, HIVpars1, mutation = TRUE, serial = TRUE) 0.02807458

##' find parameters that give specified r-valuex
find_scale <- function(target,
                       g, yini, parameters,
                       adjpar = "scale",
                       interval=c(0.9,1.3),
                       verbose=FALSE,
											 ...) {
    if (verbose) cat("target:",target,"\n")
    uu <- uniroot(function(x)
        get_rval2(x,g, yini, parameters, adjpar = adjpar,
                  verbose=verbose,...)-target,
            interval)
    return(uu$root)
}

geom_mean <- function(a){
  exp(mean(log(a)))
}
