## SHOULD BE OBSOLETE NOW
n.trial <- 1000  ## total number of sims

HIVpars.skeleton <- list(
  alphaDist=c(min=2,max=7,delta=0.5),
  scale = 1,
  scale_c = 1,
  scale_all = 1,
  Bmax = 0.317,
  Beta1 = 2.76,
  Beta3 = 0.76,
  Duration1 = 1/4,
  Duration3 = 0.75,
  K = 13938,
  a = 1.02,
  Dmax = 25.4,
  h = 0.41,
  D50 = 3058,
  Vm = 0.12,
  inisd = 0.2,
  Ini_I = 1e-4
)


HIVpars_range <- data.frame(min=c(1/10,1/15,1/5,0.01),
                            max=c(2/5,1/5,5,1),
                            row.names=c("rho_base","c_mean_base","c_u_ratio","c_e_ratio"))
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
