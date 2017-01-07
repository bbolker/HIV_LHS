library("deSolve")
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")
source("../hivModels_bd.R")

HIVpars.mean <- as.HIVvirparlist(apply(rbind(HIVpars_range, Het_range), 1, geom_mean))
HIVpars.mean <- transform(HIVpars.mean, ini_V = 3, m = 1/50)

I_mat <- matrix(NA, nrow = 1000, ncol = 8)
vir_mat <- matrix(NA, nrow = 1000, ncol = 8)
eq_vec <- rep(0, 8)
peak_mat <- matrix(NA, nrow = 8, ncol = 2)
val_vec <- rep(0, 8)

tvec <- c(1:1000)

fn <- "ev_bd.rda"

HIVpars <- HIVpars.mean
gfun1 <- gfun; gfun1_bd <- gfun_bd

fun <- function(model, i){
	yini <- calc_yini(HIVpars)
	cat("little_r...")
	little_r <- get_rval(model(HIVpars),yini, HIVpars, tvec = c(1:200))
	
	if(little_r < 0.04){
		interval <- c(1,5)
	}else{
		interval <- c(0.01,1)
	}
	cat("adj...")
	val_vec[i] <<- find_scale(0.04, model, yini, HIVpars, interval=interval, adjpar = "scale_all", tvec = c(1:200))
	HIVpars_adj <- transform(HIVpars, scale_all = val_vec[i], hmax = 0.5)
	cat("sim...\n")
	r <- rk(unlist(yini), func = model(HIVpars_adj), parms = HIVpars_adj, times = tvec, hmax = 0.5)
	
	I_mat[,i] <<- r[,(ncol(r) - 1)]
	vir_mat[,i] <<- r[,ncol(r)]
	eq_vec[i] <<- vir_mat[length(tvec),i]
	
	peak_mat[i,1] <<- which.max(vir_mat[,i])
	peak_mat[i,2] <<- max(vir_mat[,i])
	return(NULL)
}

if(file.exists(fn)){
	load(fn)
}else{
	for(i in 1:4){
		if(i == 3){
			HIVpars <- transform(HIVpars, rho_base = 5)
		}
		j <- ifelse(i %% 2 == 1, 1, 2)
		cat(2 * i - 1)
		model <- get(paste0("gfun", j))
		fun(model, 2 * i - 1)
		cat(2 * i)
		model <- get(paste0("gfun", j, "_bd"))
		fun(model, 2 * i)
		save("I_mat", "vir_mat", "eq_vec", "peak_mat", "val_vec", file = fn)
	}
}
