source("fullModel.R")
source("simFuns.R")  ## for transform.list
source("hivFuns.R")
source("Param.R")
source("hivModels.R")

HIVpars <- as.HIVvirparlist(ltab[1,])
pp <- expand(HIVpars)
yini <- calc_yini_full_ini(pp)

length(yini$S) + length(yini$I) + sum(upper.tri(yini$SS, diag = TRUE)) +
	length(yini$SI) + sum(upper.tri(yini$II, diag = TRUE))
