##
library("deSolve")
library("testthat")
source("LHSparms2.R")
source("hivFuns.R")
source("hivModels.R")

## slight asymmetry is due to enforcing row-sum=1 constraint
mm <- makeMutMat(expand(ltab_mean))
image(log10(mm))

tvec <- 0:50
fun1 <- gfun(ltab_mean,FALSE)
y1 <- unlist(calc_yini(ltab_mean))

grad1 <- fun1(t=0,yini=y1,ltab_mean)
out1 <- c(7.18942772482107e-06, -8.55762401214908e-07, 2.5370567666425e-14, 
          3.80162106809482e-10, -3.37110788855153e-07, -1.13305405627565e-05)
expect_equal(head(grad1[[1]]),out1)

res1 <- ode(y=y1,
            func=fun1,
            parms=ltab_mean,
            time=tvec)

## basic plots

## TESTS
## equivalence of models with each other?

y2 <- unlist(calc_yini(ltab2))
ltab2 <- ltab_mean
ltab2$alphaDist["delta"] <- 0.1
fun2 <- gfun(ltab2,FALSE)
res2 <- ode(y=y2,
            func=fun2,
            parms=ltab2,
            time=tvec)
