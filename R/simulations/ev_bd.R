library("deSolve")
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")
source("../hivModels_bd.R")

HIVpars.mean <- as.HIVvirparlist(apply(rbind(HIVpars_range, Het_range), 1, geom_mean))
HIVpars.mean <- transform(HIVpars.mean, ini_V = 3, m = 1/50)

I_mat_bd <- matrix(NA, nrow = 1000, ncol = 8)
vir_mat_bd <- matrix(NA, nrow = 1000, ncol = 8)
eq_vec_bd <- rep(0, 8)
peak_mat_bd <- matrix(NA, nrow = 8, ncol = 2)
val_vec_bd <- rep(0, 8)

tvec <- c(1:1000)

fn <- "ev_bd.rda"

HIVpars <- HIVpars.mean
gfun1 <- gfun; gfun1_bd <- gfun_bd

fun <- function(model, i){
	if(i %% 2 == 0){
		calc_y <- calc_yini_bd
	}else{
		calc_y <- calc_yini
	}
	
	yini <- calc_y(HIVpars)
	cat("little_r...")
	little_r <- get_rval(model(HIVpars),yini, HIVpars, tvec = c(1:200))
	
	if(little_r < 0.04){
		interval <- c(1,5)
	}else{
		interval <- c(0.01,1)
	}
	cat("adj...")
	val_vec_bd[i] <<- find_scale(0.04, model, yini, HIVpars, interval=interval, adjpar = "scale_all")
	HIVpars_adj <- transform(HIVpars, scale_all = val_vec_bd[i])
	cat("sim...\n")
	r <- rk(unlist(yini), func = model(HIVpars_adj), parms = HIVpars_adj, times = tvec, hmax = 0.1)
	
	I_mat_bd[,i] <<- r[,(ncol(r) - 1)]
	vir_mat_bd[,i] <<- r[,ncol(r)]
	eq_vec_bd[i] <<- vir_mat_bd[length(tvec),i]
	
	peak_mat_bd[i,1] <<- which.max(vir_mat_bd[,i])
	peak_mat_bd[i,2] <<- max(vir_mat_bd[,i])
	return(NULL)
}

if(file.exists(fn)){
	load(fn)
}else{
	for(i in 1:4){
		if(i == 3){
			HIVpars <- transform(HIVpars, rho_base = 5, c_u_ratio = 0)
		}
		j <- ifelse(i %% 2 == 1, 1, 2)
		cat(2 * i - 1)
		model <- get(paste0("gfun", j))
		fun(model, 2 * i - 1)
		cat(2 * i)
		model <- get(paste0("gfun", j, "_bd"))
		fun(model, 2 * i)
		save("I_mat_bd", "vir_mat_bd", "eq_vec_bd", "peak_mat_bd", "val_vec_bd", file = fn)
	}
}
