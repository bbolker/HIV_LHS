n.trial <- 1000  ## total number of sims

##Copied most descriptions from Shirreff's supplementary file. Marked as (Shirreff).
##I got the notations from Herbeck et al when I didn't know how to replicate Shirreff's code. Should we change it to Shirreff's notation?

HIVpars.skeleton <- list(
  alphaDist=c(min=2,max=7,delta=0.25),  ##distribution of strains
  scale = 1,                            ##scales transmission rate only
  scale_c = 1,                          ##scales contact rate (partnership formation/dissolution rate) only
  scale_all = 1,                        ##scales both transmission rate and contact rate
  Bmax = 0.317,                         ##Maximum rate of transmission during asymptomatic stage (Shirreff)
  K = 13938,                            ##SPVL at which infectiousness is half maximum (Shirreff)
  a = 1.02,                             ##Hill coefficient of the transmission rate
  Dmax = 25.4,                          ##Maximum duration of asymptomatic stage (Shirreff)
  h = 0.41,                             ##Hill coefficient of the duration
  D50 = 3058,                           ##SPVL at which duration of asymptomatic infection is half maximum (Shirreff)
  Vm = 0.12,                            ##Mutational standard deviation of log10SPVL (Shirreff)
  inisd = 0.2,                          ##Standard deviation of the initial SPVL distribution
  Ini_I = 1e-4,                         ##Initial number of infected individuals
  ini_V = 3.5,                          ##Initial mean SPVL
  n.risk = 10                           ##number of risk groups
)


HIVpars_range <- data.frame(min=c(1.31,1.23/12,0.413,4.81/12,1/10,1/15,1/5,0.01),
                            max=c(5.09,6/12,1.28,14/12,2/5,1/5,5,1),
                            row.names=c("Beta1","Duration1","Beta3","Duration3","rho_base","c_mean_base","c_u_ratio","c_e_ratio"))
ltab <- as.data.frame(apply(HIVpars_range,1,
                            function(x) exp(seq(log(x[1]),log(x[2]),
                                                length=n.trial))))

Het_range <- data.frame(min=c(0.01, 0.103),
												max=c(100, 1.206),
												row.names=c("kappa", "mu"))
ltab.het <- as.data.frame(apply(Het_range,1,
																function(x) exp(seq(log(x[1]),log(x[2]),
																										length=n.trial))))

set.seed(101)
ltab[] <- lapply(ltab,sample)
ltab.het[] <- lapply(ltab.het,sample)

ltab <- cbind(ltab, ltab.het)

as.HIVvirparlist <- function(x) {
  res <- append(x,HIVpars.skeleton)
  class(res) <- c("list","parlist","HIVvirparlist")
  return(res)
}

HIVpars.mean <- as.HIVvirparlist(apply(HIVpars_range, 1, geom_mean))

HIVpars.shirreff <- transform(HIVpars.mean, Beta1 = 2.76, Duration1 = 0.25, Beta3 = 0.76, Duration3 = 0.75, c_mean_base = 1.25, 
                              alphaDist = c(min = 2, max = 7, delta = 0.1))
