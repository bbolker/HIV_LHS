source("hivFuns.R")
n.trial <- 1000  ## total number of sims

## baseline parameters, not including parameters that will vary
## in the LHS process
HIVpars.skeleton <- list(
  alphaDist=c(min=2,max=7,delta=0.5),
  scale = 1,
  scale_c = 1,
  scale_all = 1,
  Bmax = 0.317,
  ## Beta1 = 2.76,
  ## Beta3 = 0.76,
  ## Duration1 = 1/4,
  ## Duration3 = 0.75,
  K = 13938,
  a = 1.02,
  Dmax = 25.4,
  h = 0.41,
  D50 = 3058,
  Vm = 0.12,
  inisd = 0.2,
  Ini_I = 1e-4
)

## parameter ranges
## also see
## https://twitter.com/noamross/status/665212230581948417
HIVpars_range <-
    data.frame(min=c(1.31,1.23/12,0.413,4.81/12,1/10,1/15,1/5,0.01),
               max=c(5.09,6/12,1.28,14/12,2/5,1/5,5,1),
               row.names=c("Beta1","Duration1","Beta3","Duration3",
               "rho_base","c_mean_base","c_u_ratio","c_e_ratio"))

## full (unscrambled) sequences
ltab <- as.data.frame(apply(HIVpars_range,1,
                            function(x) exp(seq(log(x[1]),log(x[2]),
                                                length=n.trial))))

## create LHS sample
set.seed(101)
ltab[] <- lapply(ltab,sample)

ltab_mean <- as.HIVvirparlist(apply(HIVpars_range, 1, geom_mean))

