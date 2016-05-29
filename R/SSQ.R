library("lme4")
load("../simdata/combine_ev_LHS4.rda")
sL <- sum_list[["sum_mat"]]

model = sL$model
eqV = sL$eq_vir
peakV = sL$peak_vir

lmer(eqV ~ 1 + (1|model), data = sL, REML=FALSE)
lmer(peakV ~ 1 + (1|model), data = sL, REML=FALSE) 

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