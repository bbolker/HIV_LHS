calc_yini_full_ini <- function(parameters){
	with(c(parameters),{
		c = Ini_I
		I_vec = c * I_initialize(parameters)
		yini <- list(
			S  = (1 - c) * p.risk,
			I = c(outer(I_vec, p.risk)),
			SS = matrix(0, n.risk, n.risk),
			SI = matrix(0, n.risk, n.risk * n.alpha),
			II = matrix(0, n.risk * n.alpha, n.risk *  n.alpha))
		return(yini)
	})
}

full_model_rcpp <- function(t, yini, parameters){
	tmp <- unlist(rcppFull(t, yini, parameters))
	n <- length(tmp)
	grad <- tmp[-c(n-1,n)]
	
	return(list(grad,  tot_I = tmp[n-1], mean_V = tmp[n]))
}

partner_model_rcpp <- function(t, yini, parameters){
	tmp <- unlist(rcppPartner(t, yini, parameters))
	n <- length(tmp)
	grad <- tmp[-c(n-1,n)]
	
	return(list(grad,  tot_I = tmp[n-1], mean_V = tmp[n]))
}

calc_yini_full <- function(parameters){
	yini0 <- calc_yini_full_ini(parameters)
	
	p2 <- transform(parameters, scale_all = 5)
	
	r <- rk(unlist(yini0), func = partner_model_rcpp, parms = p2, times = c(1:50))
	
	r2 <- r[50,-1]
	
	with(parameters,{
		S.ind <- 1:n.risk
		I.ind <- (n.risk+1):(n.risk+n.risk*n.alpha)
		SS.ind <- (n.risk+n.risk*n.alpha+1):(n.risk+n.risk*n.alpha+n.risk*n.risk)
		SI.ind <- (n.risk+n.risk*n.alpha+n.risk*n.risk+1):(n.risk+n.risk*n.alpha+n.risk*n.risk+n.alpha*n.risk*n.risk)
		II.ind <- (n.risk+n.risk*n.alpha+n.risk*n.risk+n.alpha*n.risk*n.risk+1):
			(n.risk+n.risk*n.alpha+n.risk*n.risk+n.alpha*n.risk*n.risk+n.alpha*n.alpha*n.risk*n.risk)
		
		yini <- list(
			S = r2[S.ind],
			I = r2[I.ind],
			SS = matrix(r2[SS.ind], n.risk, n.risk),
			SI = matrix(r2[SI.ind], n.risk, n.risk * n.alpha),
			II = matrix(r2[II.ind], n.risk * n.alpha, n.risk * n.alpha)
		)
		return(yini)
	})
}