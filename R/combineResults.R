library("ggplot2"); theme_set(theme_classic())
scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)
library("gridExtra")
library("plyr")  ## for ldply()
library("dplyr") ## for full_join()
library("reshape2")  ## for melt()
library("GGally")
library("knitr")
## devtools::install_github("lionel-/ggstance")
library("ggstance")

source("../hivFuns.R")
source("../Param.R")
HIVdefs <- c(Beta1="transmission during stage 1",
             Duration1="duration of stage 1",
             Beta3="transmission during stage 3",
             Duration3="duration of stage 3",
             rho_base="rate of partnership formation",
             c_mean_base="rate of partnership dissolution",
             c_u_ratio = "relative contact rate uncoupled",
             c_e_ratio = "relative contact rate extra-couple")

## get all the batch output and combine it into a convenient (?) format
ltab <- data.frame(run=seq(nrow(ltab)),ltab)
##  combine them/rename appropriately
run_nums <- c("pairform+epc"=15,"pairform"=16,
               "instswitch"=17,"instswitch+epc"=18,"implicit"=19, "random" = 20)
## Dushoff suggests 'instant switch' for 'real serial monogamy';
##   'implicit' instead of 'fake' ...
## each entry has
##    I_mat*:
##    vir_mat*:
##    eq_vec*: equilibrium value
##    peak_mat*: peak value (time, height)
##    val_vec*: ???
for (i in run_nums)
    L <- load(paste0("ev_LHS_res",i,".rda"))
## combine
component_names <- c("I_mat","vir_mat","eq_vec","peak_mat","val_vec")
## make one big nested list of components
comp_list <- lapply(component_names,
                    function(c) {
    res <- lapply(paste0(c,run_nums),get)
    setNames(res,names(run_nums))
})
names(comp_list) <- component_names
mat_names <- c("I_mat","vir_mat")
## reduce matrices to melted data frame of mean & quantiles
##  (we don't really need to make giant pictures with all of the
##   individual lines - ribbons with mean and 95% quantiles should
##   be sufficient)
sum_mat_fun <- function(x) {
    data.frame(
        tvec=1:nrow(x),  ## ?? better time index?
        mean=rowMeans(x),
        lwr=apply(x,1,quantile,0.025),
        upr=apply(x,1,quantile,0.975))
}
sum_list <- list()
for (m in mat_names) {
    sum_list[[m]] <- ldply(comp_list[[m]],
                           sum_mat_fun, .id="model")
}
s1 <- ldply(comp_list[["eq_vec"]],
                function(x) data.frame(run=1:length(x),eq_vir=x),
                .id="model")
s2 <- ldply(comp_list[["peak_mat"]],
                function(x) {
                   x <- setNames(data.frame(x),c("peak_time","peak_vir"))
                   data.frame(run=1:nrow(x),x)
                 },
                .id="model")
sum_list[["sum_mat"]] <- full_join(s1,s2)
## unlike merge(), full_join() preserves run order

mods <- c("pairform+epc","instswitch+epc","instswitch","pairform","implicit", "random")
fixfac <- function(x) {
    mutate(x,model=factor(model,levels=mods))
}
sum_list[["vir_mat"]] <- fixfac(sum_list[["vir_mat"]])
sum_list[["I_mat"]] <-   fixfac(sum_list[["I_mat"]])

sL <- transform(sum_list[["sum_mat"]],rel_peak=peak_vir/eq_vir)
mL <- melt(sL,id.vars=c("model","run"))

m1 <- rename(melt(sum_list[["sum_mat"]],id.vars=c("model","run")),
             sumvar=variable,sumval=value)
m2 <- rename(melt(ltab,id.vars="run"),
             LHSvar=variable,LHSval=value)
mL2 <- full_join(m1,m2)
ltab2 <- ltab[,c(2,4,5,6,7)]
names(ltab2) <- c("Beta1_adj", "Beta3_adj", "rho_ad", "c_mean_adj")
s3 <- ldply(comp_list[["val_vec"]],
                function(x) {
                   x <- data.frame(run=1:length(x),ltab2 * x)
                 },
                .id="model")

m3 <- rename(melt(s3, id.vars = c("model", "run")),
             LHSvar = variable, LHSval = value)

mL3 <- full_join(m1,m3)
s4 <- ldply(comp_list[["val_vec"]],
                function(x) {
                   x <- data.frame(run=1:length(x), c_mean_base = ltab[,7], c_mean_adj = ltab[,7] * x)
                 },
                .id="model")
sL2 <- full_join(s4, sL)

save("sum_list","fixfac","mL2","mL3","sL2",file="combineResults.rda")

