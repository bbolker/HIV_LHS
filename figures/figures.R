## setup
library(knitr)
## requires ggplot 2.2.0
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                       base_family = "Times"))
library(directlabels)
library(MASS)  ## must be before dplyr (masks 'select')
library(tidyr)
library(gridExtra)
library(plyr)  ## for ldply() -- MUST be before dplyr!
library(dplyr) ## for full_join()
## library(magrittr) ## for ???
library(reshape2)  ## for melt()
library(GGally) ## need BMB version
## devtools::install_github("bbolker/GGally")
## devtools::install_github("lionel-/ggstance")
library(ggstance)

do_png <- FALSE
do_duration <- FALSE

opts_chunk$set(echo=FALSE,error=FALSE)
if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 
zero_margin <- theme(panel.spacing=grid::unit(0,"lines"))
zero_x_margin <-
    theme(panel.spacing.x=grid::unit(0, "lines"))
## use Dark2 rather than Set1 because colour #6 of Set1 is yellow (ugh/too light)
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

source("../R/hivFuns.R")
source("../R/Param.R")

## load data
load("../simdata/combineResults.rda")

## FIXME: maybe should go in hivFuns.R ?

## setup2
orig_sum_labs <- c("peak_time","peak_vir","eq_vir","rel_vir")
new_sum_labs <- c("peak time (years)",
                  "maximum mean log10 SPVL",
                  "equilibrium mean log10 SPVL",
                  "maximum:equilibrium ratio")
m_order <- c("random", "heterogeneous","pairform+epc", "pairform",
             "instswitch+epc","instswitch" , "implicit")
fixfac2 <- function(x,atop=FALSE,newlines=FALSE) {
    if (atop) new_sum_labs <-
                  gsub("(.*) (.*)","atop(\\1,\\2)",new_sum_labs)
    factor(x,
           levels=orig_sum_labs,
           labels=new_sum_labs)
}

### Figure 1

for (i in c("","2","3")) 
    load(sprintf("../simdata/ev_LHS_resS%s.rda",i))
Ls <- load("../simdata/shirreff_data.rda")
tvec <- seq(from=1,by = 0.1,length.out=nrow(vir_matS))
dimnames(vir_matS) <- list(Time=tvec,
                           run=paste0("r=",c(0.021,0.042,0.084)))

d_littleR_m <- melt(vir_matS,value.name="Mean.VL")
shirreff_df_m <- subset(res_sum_df,run==1,select=c(run,Time, Mean.VL))
shirreff_df_m$run <- "Shirreff"
df_tot <- rbind(shirreff_df_m,d_littleR_m)
tlab <- "time (years)"
g1 <- ggplot(df_tot,aes(x=Time,y=Mean.VL,colour=run,linetype=run))+
    geom_line()+labs(x=tlab,
            y=expression(population~mean~set-point~viral~load~(log[10])))+
		scale_x_continuous(expand = c(0.005, 0), breaks = seq(0, 600, by = 100), 
											 labels = c("0", "", "200", "", "400", "", "")) +
    theme(legend.position=c(0.63,0.45),
    			legend.justification = c(0, 1),
    			panel.grid = element_blank())

## Compute difference between 

df_tot %>%
    subset(run %in% c("Shirreff", "r=0.042")) %>%
    group_by(run) %>%
    summarise(Mean.VL=max(Mean.VL)) %>%
    select(Mean.VL) %>%
    unlist -> TODO.values

## ratio of peak SPVL, Shirreff vs. single-stage model with r=0.042
(shirref_vs_single_ratio <- TODO.values[2]/TODO.values[1])

dimnames(vir_matS2) <- list(Time=tvec,
               run=c("10^{-3}", "10^{-4}", "10^{-5}"))
g2 <- g1 %+% melt(vir_matS2,id.var="Time",value.name="Mean.VL") +
	labs(y = "")
###
dimnames(vir_matS3) <- list(Time=tvec,
                            run=c("2.5", "3", "3.5"))

vir_matM3 <- melt(vir_matS3,id.var="Time",value.name="Mean.VL")
vir_matM3$run = factor(vir_matM3$run, levels = c("2.5", "3", "3.5"))

g3 <- g1 %+% vir_matM3 +
	labs(y = "")
limits <- c(2.5,4.7)
breaks <- seq(limits[1],limits[2], by = 0.5)
afun <- function(label,size=8) {
    annotate(geom="text",label=label,x=10,y=Inf,
             ## http://stackoverflow.com/questions/20083700/how-to-have-annotated-text-style-to-inherit-from-theme-set-options
             family= theme_get()$text[["family"]],
             size=size,
             vjust=1.5,hjust=-0.5)
             ## vjust=0.98,hjust=0.02)
}
g1.y <- g1 + scale_y_continuous(limits = limits, breaks = breaks)+
    scale_colour_discrete(name=bquote(run))+afun("a")+
    scale_linetype_discrete(name=bquote(run))

g2.labs <- list(bquote(10^{-3}), bquote(10^{-4}), bquote(10^{-5}))

g2.y <- g2 + scale_y_continuous(limits = limits, breaks = breaks)+
    scale_colour_discrete(name=bquote(I*"(0)"), labels = g2.labs)+afun("b")+
    scale_linetype_discrete(name=bquote(I*"(0)"),labels=g2.labs)

g3.y <- g3 + scale_y_continuous(limits = limits, breaks = breaks)+
    scale_colour_discrete(name=bquote(alpha*"(0)"))+afun("c")+
    scale_linetype_discrete(name=bquote(alpha*"(0)"))

## r fig1,fig.height=3.5,fig.width=10, echo = FALSE, dpi = 600}
fig1 <- arrangeGrob(g1.y, g2.y, g3.y, nrow=1)
ggsave("fig1.pdf", fig1, height = 3.5, width = 10)

if (do_png) {
    png(file="fig1.png",height=3.5*600,width=10*600)
    print(grid.arrange(g1.y, g2.y, g3.y, nrow=1))
    dev.off()
}

fig_objects <- c("g1.y","g2.y","g3.y")

### Figure 2

ss <- subset(sum_list[["vir_mat"]],tvec<750)

ss$model <- factor(ss$model,levels=m_order)

ss <- transform(ss,
                cmodel=droplevels(factor(gsub("+epc","",model,fixed=TRUE),
                                         levels=levels(model))))

## base plot for Figure 2 and related (duration, SPVL, beta trajectories)
gg_basetraj <- ggplot(ss,aes(tvec,mean,ymin=lwr,ymax=upr))+
    geom_line(aes(colour=model,linetype=model),lwd=1)+
    geom_ribbon(aes(fill = model), alpha=0.3)+
    labs(x="time (years)")+
    theme(axis.title.y = element_text(margin = margin(0,10,0,0)),
          axis.title.x = element_text(margin = margin(10,0,0,0)),
		legend.position = "none") +
	facet_wrap(~cmodel, nrow = 2) +
	zero_margin

## like top.points, but change hjust/vjust
top.points2 <- gapply.fun(transform(d[which.max(d$y), ],
                               ## *smaller* numbers move labels right/up
                                    hjust = 0.1, vjust = -1))

bottom.points2 <- gapply.fun(transform(d[which.min(d$y), ],
                               ## *smaller* numbers move labels right/up
                                    hjust = 0.1, vjust = 1.5))

## http://stackoverflow.com/questions/10547487/r-removing-facet-wrap-labels-completely-in-ggplot2

nostrips <- theme(strip.background = element_blank(),
                  strip.text.x = element_blank())

## r fig2, fig.width=6, fig.height = 4, echo = FALSE, cache = TRUE,dpi = 400}

gg_virtraj <- direct.label(gg_basetraj + nostrips +
    labs(y=expression(population~mean~set-point~viral~load~(log[10]))),
    method=list("top.points2","bumpup"))

ggsave(gg_virtraj,file="fig2.pdf",width=8,height=6)
if (do_png) ggsave(gg_virtraj,file="fig2.png",width=8,height=6,dpi=400)

fig_objects <- c(fig_objects,"gg_virtraj","bottom.points2","top.points2")

### Figure 2.1 - Transmission rate

ss_beta <- transform(ss, mean = returnBeta(mean,HIVpars.skeleton),
                      lwr = returnBeta(lwr,HIVpars.skeleton),
                      upr = returnBeta(upr,HIVpars.skeleton))

gg_betatraj <- direct.label(gg_basetraj %+% ss_beta + nostrips +
           labs(y=expression(mean~transmission~rate~during~asymptomatic~stage~(year^-1))),
             method=list("top.points2","bumpup"))

## scale_x_continuous(limits=c(0,2500))
## r fig2.1, fig.width=6, fig.height = 4, echo = FALSE, cache = TRUE,dpi = 400} 

ggsave(gg_betatraj,file="fig_S2_1.pdf",width=8,height=6)
if (do_png) ggsave(gg_betatraj,file="fig_S2_1.png",width=8,height=6,dpi=400)

## scale_x_continuous(limits=c(0,2500))

## r fig2.2, fig.width=6, fig.height = 4, echo = FALSE, cache = TRUE,dpi = 400} 

ss_dur <- transform(ss, mean = returnDur(mean,HIVpars.skeleton),
										lwr = returnDur(lwr,HIVpars.skeleton),
										upr = returnDur(upr,HIVpars.skeleton))

gg_durtraj <- direct.label(gg_basetraj %+% ss_dur + nostrips +
													 	labs(y="expected progression to AIDS (years)"),
													 method=list("bottom.points2","bumpup"))

ggsave(gg_durtraj,file="fig_S2_2.pdf",width=6,height=4)
if (do_png) ggsave(gg_durtraj,file="fig_S2_2.png",width=6,height=4,dpi=400)

### Figure 3

sL <- transform(sum_list[["sum_mat"]], rel_vir=peak_vir/eq_vir)
sL$model <- factor(sL$model, m_order)

## we get long tails because eq_vir has default value of 0.
## In other words, all the runs that ended up with error will have eq_vir of 0.
sL <- sL[-which(is.na(sL$peak_vir)),] %>% as_data_frame

mL <- na.omit(melt(sL,id.vars=c("model","run")))
## horrible hack for width of random-mixing violin (but doesn't help); subsample
## w <- which(mL$model=="random")
## mLw <- mL[-sample(w,size=length(w)*9/10,replace=FALSE),]
w <- with(mL,which(model=="random" & variable== "eq_vir"))
rval <- mean(mL$value[w])
mLw <- droplevels(mL[-w,])
mLw$variable <- factor(mLw$variable,
											 levels = c("peak_time", "peak_vir", "eq_vir", "rel_vir"),
											 labels = c("peak~time~(years)",
											 					 "maximum~mean~log[10]~SPVL",
											 					 "equilibrium~mean~log[10]~SPVL",
											 					 "peak:equilibrium~ratio"))

maxvir <- mL %>%
	group_by(model) %>%
	filter(variable == "peak_vir") %>%
	summarise(med=median(value),
                  lwr=quantile(value,0.025),
                  upr=quantile(value,0.975))

options(digits=3)
maxvir

gg_univ <- ggplot(mLw,aes(value,model,fill=model))+
	geom_violinh(width=1)+
	theme(legend.position = "none",
				panel.spacing.y = grid::unit(0.8, "lines")) +
	facet_wrap(~variable,scale="free_x", labeller = label_parsed)+
	guides(fill=guide_legend(reverse=TRUE))+
	labs(y="")+
	geom_point(data=data.frame(model="random",
                variable="equilibrium~mean~log[10]~SPVL",
                                   value=rval),
                   pch=22,size=3,show.legend=FALSE) +
  zero_x_margin

## r fig3, fig.width=8,fig.height=4.8, echo = FALSE, cache = TRUE,dpi = 600}

ggsave(gg_univ,file="fig3.pdf",width=8,height=4.8)
if (do_png) ggsave(gg_univ,file="fig3.png",width=8,height=4.8,dpi=600)

mLw_res <- mLw  %>%
  filter(model %in% c("implicit","pairform+epc","heterogeneous","random"))  %>%
  filter(variable %in% "maximum~mean~log[10]~SPVL")

## one-panel plot for talks
gg_univ_0 <- ggplot(mLw_res,aes(value,model,fill=model))+
	geom_violinh(width=1)+
	theme(legend.position = "none") +
	guides(fill=guide_legend(reverse=TRUE))+
    labs(x= expression(maximum~mean~log[10]~SPVL))
saveRDS(gg_univ_0,file="HIV_dur.rds")

fig_objects <- c(fig_objects,"gg_univ")

### Figure 3.1

sL.epi <- transform(sL, eq_t = returnBeta(eq_vir,HIVpars.skeleton),
                    peak_t = returnBeta(peak_vir,HIVpars.skeleton),
										eq_dur = returnDur(eq_vir,HIVpars.skeleton),
										peak_dur = returnDur(peak_vir,HIVpars.skeleton)) %>%
	select(model, run, eq_t, peak_t, eq_dur, peak_dur) %>%
	as_data_frame

maxtrans <- sL.epi %>%
	group_by(model) %>%
	summarise(med_t=median(peak_t),
						lwr_t=quantile(peak_t,0.025),
            upr_t=quantile(peak_t,0.975),
						med_d=median(peak_dur),
						lwr_d=quantile(peak_dur,0.025),
						upr_d=quantile(peak_dur,0.975))

mL.epi <- melt(sL.epi,id.vars=c("model","run"))
## horrible hack (but doesn't help); subsample
## w <- which(mL$model=="random")
## mLw <- mL[-sample(w,size=length(w)*9/10,replace=FALSE),]
w.epi.t <- with(mL.epi,which(model=="random" & variable=="eq_t"))
w.epi.d <- with(mL.epi,which(model=="random" & variable=="eq_dur"))
rval.epi.t <- mean(mL.epi$value[w.epi.t])
rval.epi.d <- mean(mL.epi$value[w.epi.d])
mLw.epi <- droplevels(mL.epi[-c(w.epi.t,w.epi.d),])
mLw.epi$variable <- factor(mLw.epi$variable,
                           levels = c("eq_t", "peak_t", "eq_dur", "peak_dur"),
                           labels = c("equilibrium~mean~transmission~rate~(year^-1)",
													 					 	"maximum~mean~transmission~rate~(year^-1)",
													 					 	"equilibrium~mean~progression~time~(years)",
													 					 	"minimum~mean~progression~time~(years)"))

gg_univ.epi <- ggplot(mLw.epi,aes(value,model,fill=model))+
	geom_violinh(width=1)+
	theme(legend.position = "none") +
	facet_wrap(~variable,scale="free_x", labeller = label_parsed)+
	guides(fill=guide_legend(reverse=TRUE))+
	labs(y="")+
	geom_point(data=data.frame(model="random",
                   variable=c("equilibrium~mean~transmission~rate~(year^-1)",
                   					 "equilibrium~mean~progression~time~(years)"),
                               value=c(rval.epi.t, rval.epi.d)),pch=22,size=3,show.legend=FALSE) +
	zero_x_margin

# http://stackoverflow.com/questions/19282897/how-to-add-expressions-to-labels-in-facet-wrap

## r fig3.1, fig.width=8,fig.height=4.8, echo = FALSE, cache = TRUE,dpi = 600}

ggsave(gg_univ.epi,file="fig_S2_3.pdf",width=8,height=4.8)
if (do_png) ggsave(gg_univ.epi,file="fig_S2_3.png",width=8,height=4.8,dpi=600)

fig_objects <- c(fig_objects,"gg_univ.epi")

## ```{r sumtab,as.is=TRUE,eval=FALSE}
if (FALSE) {
rq <- function(x,q,digits=3) {
    round(quantile(x,q),digits)
}
mLtab <- mL.epi %>% group_by(model,variable) %>%
    summarise(lwr=rq(value,0.025),
              med=rq(value,0.5),
              upr=rq(value,0.975)) %>%
    mutate(variable=abbreviate(variable,30)) %>%
    ungroup() %>%
    arrange(variable,model)
kable(mLtab)
}

## r sumvals,echo=FALSE,eval=FALSE
if (FALSE) {
ff <- function(model,variable,r=3) {
    mm <- mLtab[mLtab$model==model & mLtab$variable==variable,3:5]
    round(unlist(mm),r)
	}
	mindur_random <- ff("random","peak_dur",1)
	mindur_pfepc <-  ff("pairform+epc","peak_dur",1)
	mindur_implicit  <- ff("implicit","peak_dur",1)
	maxt_random <- ff("random","peak_t",2)
	maxt_pfepc <- ff("pairform+epc","peak_t",2)
	maxt_implicit  <- ff("implicit","peak_t",2)
	save(list=ls(pattern="^(maxt|mindur)"),file="../ms/maxminstats.RData")
}


## r fig4_func, echo = FALSE}
trim_gg <- function(gg,hack_spaces=TRUE,
										hack_parenthesis = TRUE) {
    n <- gg$nrow
    gg$nrow <- gg$ncol <- n-1
    v <- 1:n^2
    gg$plots <- gg$plots[v>n & v%%n!=0]
    gg$xAxisLabels <- gg$xAxisLabels[-n]
    gg$yAxisLabels <- gg$yAxisLabels[-1]
    ## ggpairs tries to parse labels, so this fails if they
    ## have spaces in them ... maybe there's a better way?
    if (hack_spaces) {
        gg$xAxisLabels <- gsub("_"," ",gg$xAxisLabels)
        gg$yAxisLabels <- gsub("_"," ",gg$yAxisLabels)
    }
    if (hack_parenthesis) {
    	gg$xAxisLabels <- gsub("years","\\(years\\)",gg$xAxisLabels)
    	gg$yAxisLabels <- gsub("years","\\(years\\)",gg$yAxisLabels)
    }
    
    return(gg)
}

tweak_legends_gg <- function(gg,pos=1.5) {
    n <- gg$nrow
    for (i in 1:n) {
        for (j in 1:n) {
            inner <- getPlot(gg,i,j)
            if (i==1 & j==1) {
                inner <- inner + ggplot2::theme(legend.position=c(pos,0.5)) +
                    guides(colour = guide_legend(override.aes = list(size=6)))
            } else {
                inner <- inner + ggplot2::theme(legend.position="none")
            }
            gg <- putPlot(gg,inner,i,j)
        }
    }
    return(gg)
}
tweak_colours_gg <- function(gg) {
    n <- gg$nrow
    for (i in 1:n) {
        for (j in 1:n) {
            inner <- getPlot(gg,i,j)
            inner <- inner+scale_colour_brewer(palette = 'Dark2')+
                scale_shape_manual(values=1:7)
            gg <- putPlot(gg,inner,i,j)
        }
    }
    return(gg)
}

### Figure 4

sL2 <- sL[,c("model", "run", "eq_vir", "peak_time", "peak_vir")]

set.seed(101)
sL2 <- sL2 %>%
	group_by(model) %>%
	sample_n(100)

new_sum_labs2 <- gsub("\\(years\\)", "years", new_sum_labs[1:3])

names(sL2)[na.omit(match(orig_sum_labs,names(sL2)))] <- gsub(" ","_",new_sum_labs2)
    ## paste0("`",new_sum_labs,"`")
ggp1 <- ggpairs(sL2,
        mapping = ggplot2::aes(color = model,pch=model),
        columns=3:5,
        legends=TRUE,
        lower = list(continuous = wrap("points",alpha=0.6,size=2)), ## alpha = 0.3,size=0.5)),
        diag = list(continuous = "blankDiag"),
        upper = list(continuous = "blank"))
## http://stackoverflow.com/questions/14711550/is-there-a-way-to-change-the-color-palette-for-ggallyggpairs-using-ggplot
## have to change plot one panel at a time
ggp2 <- tweak_colours_gg(tweak_legends_gg(trim_gg(ggp1)))
## add 1-to-1 line
inner <- getPlot(ggp2,2,1)
inner <- inner+geom_abline(intercept=0,slope=1,lty=2)
ggp2 <- putPlot(ggp2,inner,2,1)
## tweak scale
for(i in 1:2){
	inner <- getPlot(ggp2,2,i)
	inner <- inner+scale_y_continuous(breaks=seq(3,5,by=1))
	ggp2 <- putPlot(ggp2,inner,2,i)	
}

## ``{r fig4,fig.width=7,fig.height=7, echo = FALSE, cache = TRUE,dpi = 600}

pdf(file="fig4.pdf",width=7,height=7)
print(ggp2,spacingProportion=0, left = 0.14)
dev.off()

if (do_png) {
    png(file="fig4.png",width=7*600,height=7*600)
    print(ggp2,spacingProportion=0)
    dev.off()
}

fig_objects <- c(fig_objects,"ggp2")

### Figure 5
remove_runs <- unique(subset(mL2, sumvar == "peak_time" & is.na(sumval))[,"run"])
mL2 <- subset(mL2, !(model == "heterogeneous" & run %in% remove_runs))

mL3 <- subset(mL2,
       !((model=="implicit" &
          LHSvar %in% c("rho_base","c_u_ratio","c_e_ratio", "kappa", "mu")) |
         (model=="random" &
          LHSvar %in% c("rho_base","c_u_ratio","c_e_ratio", "kappa", "mu")) |
         (model=="instswitch" &
          LHSvar %in% c("rho_base","c_u_ratio","c_e_ratio", "kappa", "mu")) |
         (model=="instswitch+epc" &
          LHSvar %in% c("rho_base","c_u_ratio", "kappa", "mu")) |
         (model=="pairform" &
          LHSvar %in% c("c_u_ratio","c_e_ratio", "kappa", "mu")) |
       		(model == "pairform+epc" &
       		 	LHSvar %in% c("kappa", "mu"))))

mL3 <- droplevels(mL3)
fvals <- c("Duration1","c_mean_base","c_e_ratio","rho_base",
                              "c_u_ratio",
                              "kappa", "mu")

mL3 <- subset(mL3,LHSvar %in% fvals)

## reorder factors
mL3$LHSvar <- factor(mL3$LHSvar,
                     levels=fvals,
                     labels = c("D[P]",
                                "c", "c[e]/c[w]", "rho", "c[u]/c[w]",
                                "kappa", "mu"))

## drop rows with variables we don't want anyway
mL3 <- na.omit(mL3) 

## FIXME: could also use labels argument to rename levels, labeller=label_parsed
## in facet_grid to get pretty math

mL3$model <- factor(mL3$model, levels=m_order)

mL3 <- subset(mL3,sumvar %in% orig_sum_labs)
mL3$sumvar <- fixfac2(mL3$sumvar)

L <- function(labels,multi_line=TRUE) {
    r <- if (all(grepl(" ",labels[[1]]))) {
             list(as.character(gsub(" ","\n",labels[[1]])))
         } else {
             label_parsed(labels,multi_line=multi_line)
         }
    return(r)
}

brkfun <- function(x) {
    nx <- 1
    n <- 2
    brks <- axisTicks(log(x), log=TRUE,  nint = n)
    ## browser()
    brks <- brks[brks>x[1] & brks<x[2]]
    nx <- length(brks)
    ## cat("brk1",brks,"\n")
    while (nx<2 && n<10) {
        n <- n+1
        ## cat(n,"\n")
        brks <- axisTicks(log10(x), log=TRUE,  nint = n)
        brks <- brks[brks>x[1] & brks<x[2]]
        ## cat("brk2",brks,"\n")
        nx <- length(brks)
    }
    if (n==10) browser()
    r <- brks[c(1,length(brks))]
    ## hacks.
    ## would like to customize x-axis labels for each column,
    ##  but that's not so easy
    ## So instead is where we try to adjust the last few axis ticks ...
    ## r needs to be rounded if we want to compare it.. probably long values causing problems...
    r <- round(r, digit = 2)
    if (all(r==c(0.01,1))) r <- c(0.02,0.5)
    if (all(r==c(0.2,5))) r <- c(0.5,5)
    if (all(r==c(0.01,100))) r <- c(0.05, 50)
    if (all(r==c(0.1,1))) r <- c(0.5, 1)
    return(r)
}

## r fig5, echo = FALSE, fig.width = 10, cache = TRUE,dpi = 600}

ggsens <- ggplot(mL3,aes(LHSval,sumval,colour=model))+
    geom_point(pch=".",alpha=0.5)+
    facet_grid(sumvar~LHSvar,scales="free",
               labeller = L)+
    geom_smooth(se=FALSE)+labs(x="",y="")+
    ## leave a little extra room?
    scale_x_log10(expand=c(0,0.08), breaks=brkfun)+
		scale_y_continuous(expand=c(0,0)) +
    zero_margin

ggsave(ggsens, file="fig5.pdf",width=10,height=5)
if (do_png) ggsave(ggsens, file="fig5.png",width=10,height=5,dpi=600)

fig_objects <- c(fig_objects,"ggsens")

save(list=c("fig_objects",fig_objects),file="HIVLHS_figures.RData")
