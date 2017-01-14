library("Rcpp")
library("deSolve")
sourceCpp("../fullModel.cpp")
source("../fullModel.R")
source("../simFuns.R")  ## for transform.list
source("../hivFuns.R")
source("../Param.R")
source("../hivModels.R")

HIVpars.mean <- as.HIVvirparlist(apply(rbind(HIVpars_range, Het_range), 1, geom_mean))
HIVpars.mean <- transform(HIVpars.mean, ini_V = 3)

param_list <- setNames(rep(list(HIVpars.mean), 4),
                       c("pairform+epc","pairform",
                         "instswitch+epc","instswitch"))

## set large pair-formation rate for 'instswitch' models
param_list[["instswitch+epc"]]$rho_base <-
    param_list[["instswitch"]]$rho_base <- 4

## set no extrapair/no-uncoupled contact for 'noepc' models
extra_un_params <- c("c_u_ratio", "c_e_ratio")
param_list[["pairform"]][extra_un_params] <- 
	param_list[["instswitch"]][extra_un_params] <- 
	c(0, 0)

param_list[["instswitch+epc"]]$c_u_ratio <- 0

tvec <- 1:1000

I_matFull <- matrix(NA, nrow = length(tvec), ncol = 8)
vir_matFull <- matrix(NA, nrow = length(tvec), ncol = 8)
eq_vecFull <- rep(0, 8)
peak_matFull <- matrix(NA, nrow = 8, ncol = 2)
val_vecFull <- rep(0, 8)

fn <- "ev_hetero.rda"

fun <- function(pars, model, i) {
    HIVpars <- pars
    ## pair <- (i<=4) ?
    if (i > 4) {
        pair = FALSE
    } else {
        pair = TRUE
    }
    pp <- expand(HIVpars, pair = pair)
    yini <- calc_yini_full(pp)
    cat("finding r...")
    little_r <- get_rval(model, yini, pp, hmax = 0.5, tvec = c(1:200))
    
    if (little_r < 0.04) {
        interval <- c(1, 10)
    } else {
        interval <- c(0.01, 1)
    }
    cat("adj...")
    val_vecFull[i] <<- find_scale(0.04, model, yini, HIVpars, interval=interval,
                                  adjpar = "scale_all", full = TRUE, hmax = 0.5, tvec = 1:400, verbose = TRUE)
    HIVpars_adj <- transform(HIVpars, scale_all = val_vecFull[i])
    pp_adj <- expand(HIVpars_adj)
    cat("sim...\n")
    r <- rk(unlist(yini), func = model, parms = pp_adj, times = tvec, hmax = 0.1)
    
    I_matFull[,i] <<- r[,(ncol(r) - 1)]
    vir_matFull[,i] <<- r[,ncol(r)]
    eq_vecFull[i] <<- vir_matFull[length(tvec),i]
    
    peak_matFull[i,1] <<- which.max(vir_matFull[,i])
    peak_matFull[i,2] <<- max(vir_matFull[,i])
    return(NULL)
}

if (file.exists(fn)) {
    load(fn)
} else {
    for (i in 1:4) {
        model <- ifelse(i == 1 | i == 3, full_model_rcpp, full_model_rcpp_noepc)
        ## odd values are for no-heterogeneity
        cat(2 * i - 1)
        base_pars <- param_list[[i]]
        fun(base_pars, model, 2 * i - 1)
        cat(2 * i)
        nohetero_pars <- transform(base_pars, n.risk = 1)
        fun(nohetero_pars, model, 2 * i)
        save("I_matFull", "vir_matFull", "eq_vecFull", "peak_matFull", "val_vecFull", file = fn)
    }
}
