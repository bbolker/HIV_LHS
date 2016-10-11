simdatapath <- "../simdata"
source("hivFuns.R")  ## must be loaded before LHSparms3
source("LHSparms3.R")
ltab <- data.frame(run=seq(nrow(ltab)),ltab)
##  combine them/rename appropriately
run_nums <- c("pairform+epc"=5,"pairform"=6,
              "instswitch"=7,"instswitch+epc"=8,"implicit"=9)
####              "heterogeneous"=10)
## each entry has
##    I_mat*:
##    vir_mat*:
##    eq_vec*: equilibrium value
##    peak_mat*: peak value (time, height)
##    val_vec*: ???
for (i in run_nums)
    L <- load(file.path(simdatapath,paste0("ev_LHS_res",i,".rda")))
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
```

```{r fix_order}
mods <- c("pairform+epc","instswitch+epc","instswitch","pairform","implicit")
fixfac <- function(x) {
    mutate(x,model=factor(model,levels=mods))
}
sum_list[["vir_mat"]] <- fixfac(sum_list[["vir_mat"]])
sum_list[["I_mat"]] <-   fixfac(sum_list[["I_mat"]])
