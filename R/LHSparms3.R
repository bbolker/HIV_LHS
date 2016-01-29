n.trial <- 1000  ## total number of sims

## FIXME: describe all of these!
HIVpars.skeleton <- list(
  alphaDist=c(min=2,max=7,delta=0.25),
  scale = 1,
  scale_c = 1,
  scale_all = 1,
  Bmax = 0.317,
  K = 13938,
  a = 1.02,
  Dmax = 25.4,
  h = 0.41,
  D50 = 3058,
  Vm = 0.12,
  inisd = 0.2,
  Ini_I = 1e-4,
  ini_V = 3.5
)


## FIXME: document all of these!
HIVpars_range <- data.frame(min=c(1.31,1.23/12,0.413,
                                  4.81/12,1/10,1/15,1/5,0.01),
                            max=c(5.09,6/12,1.28,14/12,2/5,1/5,5,1),
                            row.names=c("Beta1","Duration1",
                                        "Beta3","Duration3",
                                        "rho_base","c_mean_base",
                                        "c_u_ratio","c_e_ratio"))


ltab <- as.data.frame(apply(HIVpars_range,1,
                            function(x) exp(seq(log(x[1]),log(x[2]),
                                                length=n.trial))))

set.seed(101)
ltab[] <- lapply(ltab,sample)

as.HIVvirparlist <- function(x) {
  res <- append(x,HIVpars.skeleton)
  class(res) <- c("list","parlist","HIVvirparlist")
  return(res)
}

ltab_mean <- as.HIVvirparlist(apply(HIVpars_range, 1, geom_mean))

