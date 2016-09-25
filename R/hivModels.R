#extra couple
gfun <- function(parameters, experimental=FALSE) {
	px <- expand(parameters)
  g <- function(t,yini,parameters) {
      with(as.list(c(yini,px)), 
      { 
          SI <- yini[3:(2+n.alpha)]
          I <- yini[(3+n.alpha):(2+2*n.alpha)]
          II <- matrix(yini[(3+2*n.alpha):(2+(2+n.alpha)*n.alpha)],
                       n.alpha, n.alpha)
  
          II_adj <- II
          diag(II_adj) <- 2 * diag(II_adj) 
  
          prop <- (c_u_ratio * I + c_e_ratio * (SI + colSums(II_adj)))/
              (c_u_ratio * (S + sum(I)) +
               c_e_ratio * (2 * SS + 2 * sum(SI) + sum(II_adj)))
  
          dS <- - rho * S * S/(S+sum(I)) +
              2 * c_mean * SS - rho * S *sum( I)/(S+sum(I)) +
              c_mean * sum(SI) + 2 * sum(lam * SI) + sum(lam * I) +
              sum(lammat_dis * II) - sum(prop * c_u) * S
          dSS <- 0.5 * rho * S * S/(S+sum(I)) - c_mean * SS -
              2 * sum(prop * c_e) * SS
          dSI <- rho * S * I/(S+sum(I)) - c_mean * SI -
              lam * SI - Beta * SI +
              (prop * 2 * c_e * SS) %*% p - sum(c_e * prop) * SI
  
          infrate <- (Beta * SI) * p
          infrate_adj <- infrate + t(infrate)
          diag(infrate_adj) <- diag(infrate_adj)/2
  
          infrate_e <- outer(SI, c((c_e * prop) %*% p))
          infrate_e_adj <- infrate_e + t(infrate_e)
          infrate_e_save <- infrate_e_adj
          diag(infrate_e_adj) <- diag(infrate_e_adj)/2
  
          frate.I <- rho * outer(I,I,"*")/(S+sum(I))
          diag(frate.I) <- diag(frate.I)/2
  
          dI <- - rho * S * I/(S+sum(I)) +
              c_mean * SI - rho * I * sum(I)/(S+sum(I))  +
              colSums(c_mean * II_adj) - lam * I +
              colSums(lammat_dis * II) + (prop * c_u * S) %*% p
          dII <- frate.I - c_mean * II - lammat_adj * II +
              infrate_adj + infrate_e_adj
  
          tot <- S + 2* sum(SI) + 2*SS + sum(I) + sum(II_adj)
  
          tot_I <- sum(SI) + sum(I) + sum(II_adj)
  
          tot_V <- sum(SI * alpha) + sum(I * alpha) +
              sum(sweep(II_adj, 2, alpha, "*"))
          list(c(dS, dSS, dSI, dI, dII), tot_I = tot_I, mean_V = tot_V/tot_I)
      })
  }

  ## tweak for performance
	gX <- function(t,yini,parameters) {
      with(as.list(c(yini,px)),
      { 
          SI <- yini[3:(2+n.alpha)]
          I <- yini[(3+n.alpha):(2+2*n.alpha)]
          II <- matrix(yini[(3+2*n.alpha):(2+(2+n.alpha)*n.alpha)],
                       n.alpha, n.alpha)
  
          II_adj <- II
          diag(II_adj) <- 2 * diag(II_adj) 

          ## precompute sums for efficiency
          sumI <- sum(I)
          sumSI <- sum(SI)
          
          prop <- (c_u_ratio * I + c_e_ratio * (SI + colSums(II_adj)))/
              (c_u_ratio * (S + sumI) +
               c_e_ratio * (2 * SS + 2 * sumSI + sum(II_adj)))

          sumpropce <- sum(prop*c_e)
          s <- seq_len(n.alpha)
          mind <- cbind(s,s)

          dS <- - rho * S * S/(S+sumI) +
              2 * c_mean * SS - rho * S *sumI/(S+sumI) +
              c_mean * sumSI + 2 * sum(lam * SI) + sum(lam * I) +
              sum(lammat_dis * II) - sum(prop * c_u) * S
          dSS <- 0.5 * rho * S * S/(S+sumI) - c_mean * SS -
              2 * sumpropce * SS
          dSI <- rho * S * I/(S+sumI) - c_mean * SI -
              lam * SI - Beta * SI +
              (prop * 2 * c_e * SS) %*% p - sumpropce * SI
  
          infrate <- (Beta * SI) * p
          infrate_adj <- infrate + t(infrate)
          infrate_adj[mind] <- infrate_adj[mind]/2
  
          infrate_e <- outer(SI, c((c_e * prop) %*% p))
          infrate_e_adj <- infrate_e + t(infrate_e)
          infrate_e_save <- infrate_e_adj
          infrate_e_adj[mind] <- infrate_e_adj[mind]/2
  
          frate.I <- rho * outer(I,I,"*")/(S+sumI)
          frate.I[mind] <- frate.I[mind]/2
  
          dI <- -rho * S*I/(S+sumI) +
              c_mean * SI - rho*I*sumI/(S+sumI)  +
              colSums(c_mean*II_adj) - lam*I +
              colSums(lammat_dis*II) + (prop * c_u * S) %*% p
          dII <- frate.I - c_mean*II - lammat_adj*II +
              infrate_adj + infrate_e_adj
  
          tot <- S + 2*sumSI + 2*SS + sumI + sum(II_adj)
          tot_I <- sumSI + sumI + sum(II_adj)
          tot_V <- sum(SI * alpha) + sum(I * alpha) +
              sum(sweep(II_adj, 2, alpha, "*"))
          
          list(c(dS, dSS, dSI, dI, dII), tot_I = tot_I, mean_V = tot_V/tot_I)
      })
  } 

  if (experimental) return(gX) else return(g)
}

#serial
gfun2 <- function(parameters) {
	px <- expand(parameters)
  g <- function(t,yini,parameters) {
    with(as.list(c(yini,px)), 
{      
  
  SI = yini[3:(2+n.alpha)]
  I = yini[(3+n.alpha):(2+2*n.alpha)]
  II = matrix(yini[(3+2*n.alpha):(2+(2+n.alpha)*n.alpha)], n.alpha, n.alpha)
  
  dS = - rho * S + 2 * c_mean * SS + c_mean * sum(SI) + 2 * sum(lam * SI) + sum(lam * I) + sum(lammat_dis * II)
  dSS = 0.5 * rho * S * S/(S+sum(I)) - c_mean * SS
  dSI = rho * S * I/(S+sum(I)) - c_mean * SI - lam * SI - Beta * SI
  
  infrate = (Beta * SI) * p
  infrate_adj = infrate + t(infrate)
  diag(infrate_adj) = diag(infrate_adj)/2
  
  frate.I = rho * outer(I,I,"*")/(S+sum(I))
  diag(frate.I) = diag(frate.I)/2
  
  II_adj = II
  diag(II_adj) = 2 * diag(II_adj) 
  
  dI = - rho * I + c_mean * SI + colSums(c_mean * II_adj) - lam * I + colSums(lammat_dis * II)
  dII = frate.I - c_mean * II - lammat_adj * II + infrate_adj
  
  tot = S + 2* sum(SI) + 2*SS + sum(I) + sum(II_adj)
  
  tot_I = sum(SI) + sum(I) + sum(II_adj)
  
  tot_V = sum(SI * alpha) + sum(I * alpha) + sum(sweep(II_adj, 2, alpha, "*"))
  
  list(c(dS, dSS, dSI, dI, dII), tot_I = tot_I, mean_V = tot_V/tot_I)
})
  } 
return(g)
}

#serial + inst
gfun3 <- function(parameters) {
	px <- expand(parameters)
  g <- function(t,yini,parameters) {
    with(as.list(c(yini,px)), 
{      
  
  SI = yini[2:(1+n.alpha)]
  II = matrix(yini[(2+n.alpha):(1+(1+n.alpha)*n.alpha)], n.alpha, n.alpha)
  
  
  II_adj = II
  diag(II_adj) = 2 * diag(II_adj)
  
  S = 2 * c_mean * SS + c_mean * sum(SI) + 2 * sum(lam * SI) + sum(lammat_dis * II)
  
  II_adj = II
  diag(II_adj) = 2 * diag(II_adj)
  
  I = c_mean * SI + colSums(c_mean * II_adj) + colSums(lammat_dis * II)
  
  dSS = - c_mean * SS + 0.5 * S * S /(S + sum(I))
  
  dSI = - c_mean * SI - Beta * SI - lam * SI + S * I/(S+sum(I))
  
  infrate = (Beta * SI) * p
  infrate_adj = infrate + t(infrate)
  diag(infrate_adj) = diag(infrate_adj)/2
  
  f.rate = outer(I, I, "*")/(S + sum(I))
  diag(f.rate) = diag(f.rate)/2
  
  dII = - c_mean * II + infrate_adj - lammat_adj * II + f.rate
  
  
  tot = 2* sum(SI) + 2*SS + sum(II_adj)
  
  tot_I = sum(SI) + sum(II_adj)
  
  tot_V = sum(SI * alpha) + sum(sweep(II_adj, 2, alpha, "*"))
  
  list(c(dSS, dSI, dII), tot_I = tot_I, mean_V = tot_V/tot_I)
})
  } 
return(g)
}


#serial + inst + extra
gfun4 <- function(parameters) {
	px <- expand(parameters)
  g <- function(t,yini,parameters) {
    with(as.list(c(yini,px)), 
{      
  
  SI = yini[2:(1+n.alpha)]
  II = matrix(yini[(2+n.alpha):(1+(1+n.alpha)*n.alpha)], n.alpha, n.alpha)
  
  
  II_adj = II
  diag(II_adj) = 2 * diag(II_adj)
  prop = (c_e_ratio * (SI + colSums(II_adj)))/(c_e_ratio * (2 * SS + 2 * sum(SI) + sum(II_adj)))
  S = 2 * c_mean * SS + c_mean * sum(SI) + 2 * sum(lam * SI) + sum(lammat_dis * II)
  
  II_adj = II
  diag(II_adj) = 2 * diag(II_adj)
  
  I = c_mean * SI + colSums(c_mean * II_adj) + colSums(lammat_dis * II)
  
  dSS = - c_mean * SS + 0.5 * S * S /(S + sum(I)) - 2 * sum(prop * c_e) * SS
  
  dSI = - c_mean * SI - Beta * SI - lam * SI + S * I/(S+sum(I)) + (prop * 2 * c_e * SS) %*% p - sum(c_e * prop) * SI
  
  infrate_e = outer(SI, c((c_e * prop) %*% p))
  infrate_e_adj = infrate_e + t(infrate_e)
  infrate_e_save = infrate_e_adj
  diag(infrate_e_adj) = diag(infrate_e_adj)/2
  
  infrate = (Beta * SI) * p
  infrate_adj = infrate + t(infrate)
  diag(infrate_adj) = diag(infrate_adj)/2
  
  f.rate = outer(I, I, "*")/(S + sum(I))
  diag(f.rate) = diag(f.rate)/2
  
  dII = - c_mean * II + infrate_adj - lammat_adj * II + f.rate + infrate_e_adj
  
  
  tot = 2* sum(SI) + 2*SS + sum(II_adj)
  
  tot_I = sum(SI) + sum(II_adj)
  
  tot_V = sum(SI * alpha) + sum(sweep(II_adj, 2, alpha, "*"))
  
  list(c(dSS, dSI, dII), tot_I = tot_I, mean_V = tot_V/tot_I)
})
  } 
return(g)
}

#shirreff
gfun5 <- function(parameters) {
	px <- expand(parameters)
  g <- function(t,yini,parameters) {
    with(as.list(c(yini,px)), 
{      
  
  I = yini[-1]
  
  Beta_adj = c_mean * Beta/(c_mean + Beta + lam)
  
  dS = sum(lam * I) - sum(Beta_adj * I * S)

  dI = (Beta_adj * I * S) %*% p - lam * I
  
  tot_I = sum(I)
  
  tot_V = sum(alpha * I)
  
  list(c(dS, dI), tot_I = tot_I, mean_V = tot_V/tot_I)
})
  } 
return(g)
}


#random mixing
gfun_random <- function(parameters) {
	px <- expand(parameters)
  g <- function(t,yini,parameters) {
    with(as.list(c(yini,px)), 
{      
  
  I = yini[-1]
 
  dS = sum(lam * I) - sum(c_mean * Beta * I * S)
  
  dI = (c_mean * Beta * I * S) %*% p - lam * I
  
  tot_I = sum(I)
  
  tot_V = sum(alpha * I)
  
  list(c(dS, dI), tot_I = tot_I, mean_V = tot_V/tot_I)
})
  } 
return(g)
}


I_initialize <- function(parameters) {
  with(expand(parameters),
  {
    I = dnorm(alpha, mean = ini_V, sd = inisd)/sum(dnorm(alpha, mean = ini_V, sd = inisd))
    return(I)
  })
}

calc_yini <- function(parameters){
  with(c(expand(parameters)),{
    x = (c_mean)/(c_mean + rho)
    N = (1-x)/2
    c = Ini_I
    I_vec = c * I_initialize(parameters)
    I_mat = 2 * outer(I_vec, I_vec)
    diag(I_mat) = diag(I_mat)/2
    yini <- list(
      S  = (1 - c) * x,
      SS = (1 - c)^2 * N,
      SI = 2 * (1-c) * N * I_vec,
      I = x * I_vec,
      II = N * I_mat)
    return(yini)
  })
}

calc_yini2 <- function(parameters){
  with(c(expand(parameters)),{
    N = 0.5
    c = Ini_I
    I_vec = c * I_initialize(parameters)
    I_mat = 2 * outer(I_vec, I_vec)
    diag(I_mat) = diag(I_mat)/2
    yini <- list(
      SS = (1 - c)^2 * N,
      SI = 2 * (1-c) * N * I_vec,
      II = N * I_mat)
    return(yini)
  })
}

calc_yini3 <- function(parameters){
  with(c(expand(parameters)),{
    yini <- list(
      S = 1 - Ini_I,
      I = Ini_I * I_initialize(parameters)
    )
    return(yini)
  })
}

tvec <- seq(0,4000,by=1)
