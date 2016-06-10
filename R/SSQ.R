library("lme4")
load("../simdata/combine_ev_LHS4.rda")
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

##why is this different from the lmer results?

## do some calculations on the scales of transmission prob/duration as
## well as on the log10 SPVL scale, for example.
##
##    calculate the peak transmission probability and the 95% quantiles
##    of peak transmission probability for the pairform+epc model
##
##    ditto for the implicit model

##    ditto for duration
