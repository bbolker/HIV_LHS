library("lme4")
source("hivFuns.R")
source("Param.R")
load("../simdata/combine_ev_LHS4.rda")
source("hivFuns.R")
source("Param.R")

sL <- sum_list[["sum_mat"]]

model = sL$model
eqV = sL$eq_vir
peakV = sL$peak_vir

m_peakV <- lmer(peakV ~ 1 + (1|model), data = sL, REML=FALSE)

v_among <- c(VarCorr(m_peakV)[[1]])  ## c() dumps extra information
v_within <- sigma(m_peakV)^2
v_among/(v_among+v_within)


library(nlme)
lme(eqV~1,random=~1|model,data=sL,method="ML")

findSSQ <- function(var){
	wSSQ = sum((var-ave(var, model))^2) ##within
	group_mean = tapply(var,list(model),mean)
	bSSQ = sum(((group_mean)-mean(var))^2)*1000 ##between
	
	cat("within model standard deviation:", sqrt(wSSQ/6000),"\n")
	cat("between model standard deviation:", sqrt(bSSQ/6000),"\n")
}

findSSQ(eqV)
findSSQ(peakV)

## do some calculations on the scales of transmission prob/duration as
## well as on the log10 SPVL scale, for example.
##
##    calculate the peak transmission probability and the 95% quantiles
##    of peak transmission probability for the pairform+epc model
##
##    ditto for the implicit model

##    ditto for duration


sL <- transform(sum_list[["sum_mat"]],rel_peak=peak_vir/eq_vir)
sL$model <- factor(sL$model, c("random","pairform+epc", "pairform", "instswitch+epc","instswitch" , "implicit"))

returnBeta <- function(lSpvl, parameters){
	with(as.list(c(parameters)),{
		v = 10^lSpvl
		return(hill(K,Bmax,v,a))
	})
}

returnDur <- function(lSpvl, parameters){
	with(as.list(c(parameters)),{
		v = 10^lSpvl
		return(hill(v,Dmax,D50,h))
	})
}

sum_fun <- function(x) {
	data.frame(
		var = names(x),
		mean=colMeans(x),
		lwr=apply(x,2,quantile,0.025),
		upr=apply(x,2,quantile,0.975))
}

sL.epi <- transform(sL, eq_t = returnBeta(eq_vir,HIVpars.skeleton),
										peak_t = returnBeta(peak_vir,HIVpars.skeleton),
										eq_dur = returnDur(eq_vir,HIVpars.skeleton),
										peak_dur = returnDur(peak_vir,HIVpars.skeleton))
sL.epi <- sL.epi[,-c(3:6)]

sL.epi

model_name <- c("pairform+epc", "pairform", "instswitch", "instswitch+epc", "implicit", "random")

epiList <- list()

for(i in 1:6){
	epiList[[model_name[i]]] = sL.epi[c((1000*(i-1)+1):(1000*i)),3:6]
}

library("plyr")

ldply(epiList, sum_fun)
