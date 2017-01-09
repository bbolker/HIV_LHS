gfun_bd <- function(parameters, experimental=FALSE) {
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
				 		c_mean * sum(SI) + sum(lam * SI) - sum(prop * c_u) * S +
				 		m * (1 - S) + 2 * m * SS + m * sum(SI)
				 	dSS <- 0.5 * rho * S * S/(S+sum(I)) - c_mean * SS -
				 		2 * sum(prop * c_e) * SS - 2 * m * SS
				 	dSI <- rho * S * I/(S+sum(I)) - c_mean * SI -
				 		lam * SI - Beta * SI +
				 		(prop * 2 * c_e * SS) %*% p - sum(c_e * prop) * SI - 2 * m * SI
				 	
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
				 		colSums(lammat_dis * II) + (prop * c_u * S) %*% p -
				 		m * I + m * SI + colSums(m * II_adj)
				 	dII <- frate.I - c_mean * II - lammat_adj * II +
				 		infrate_adj + infrate_e_adj - 2 * m * II
				 	
				 	tot <- S + 2* sum(SI) + 2*SS + sum(I) + sum(II_adj)
				 	
				 	tot_I <- sum(SI) + sum(I) + sum(II_adj)
				 	
				 	tot_V <- sum(SI * alpha) + sum(I * alpha) +
				 		sum(sweep(II_adj, 2, alpha, "*"))
				 	list(c(dS, dSS, dSI, dI, dII), tot = tot, I = tot_I/tot, mean_V = tot_V/tot_I)
				 })
	}
	return(g)
}

#serial
gfun2_bd <- function(parameters) {
	px <- expand(parameters)
	g <- function(t,yini,parameters) {
		with(as.list(c(yini,px)), 
				 {      
				 	
				 	SI = yini[3:(2+n.alpha)]
				 	I = yini[(3+n.alpha):(2+2*n.alpha)]
				 	II = matrix(yini[(3+2*n.alpha):(2+(2+n.alpha)*n.alpha)], n.alpha, n.alpha)
				 	
				 	II_adj = II
				 	diag(II_adj) = 2 * diag(II_adj) 
				 	
				 	dS = - rho * S + 2 * c_mean * SS + c_mean * sum(SI) + sum(lam * SI) +
				 		m * (1 - S) + 2 * m * SS + m * sum(SI)
				 	dSS = 0.5 * rho * S * S/(S+sum(I)) - c_mean * SS - 2 * m * SS
				 	dSI = rho * S * I/(S+sum(I)) - c_mean * SI - lam * SI - Beta * SI - 2 * m * SI
				 	
				 	infrate = (Beta * SI) * p
				 	infrate_adj = infrate + t(infrate)
				 	diag(infrate_adj) = diag(infrate_adj)/2
				 	
				 	frate.I = rho * outer(I,I,"*")/(S+sum(I))
				 	diag(frate.I) = diag(frate.I)/2
				 	
				 	dI = - rho * I + c_mean * SI + colSums(c_mean * II_adj) - lam * I + colSums(lammat_dis * II) - 
				 		m * I + m * SI + colSums(m * II_adj)
				 	dII = frate.I - c_mean * II - lammat_adj * II + infrate_adj - 2 * m * II
				 	
				 	tot = S + 2* sum(SI) + 2*SS + sum(I) + sum(II_adj)
				 	
				 	tot_I = sum(SI) + sum(I) + sum(II_adj)
				 	
				 	tot_V = sum(SI * alpha) + sum(I * alpha) + sum(sweep(II_adj, 2, alpha, "*"))
				 	
				 	list(c(dS, dSS, dSI, dI, dII), tot = tot, I = tot_I/tot, mean_V = tot_V/tot_I)
				 })
	} 
	return(g)
}

calc_yini_bd <- function(parameters){
	with(c(expand(parameters)),{
		x = (c_mean + 2 * m)/(c_mean + rho + 2 * m)
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
