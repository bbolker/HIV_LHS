source("fig_setup.R")

sL <- transform(sum_list[["sum_mat"]],rel_peak=peak_vir/eq_vir)
sL$model <- factor(sL$model, c("random","pairform+epc", "pairform", "instswitch+epc","instswitch" , "implicit"))
mL <- melt(sL,id.vars=c("model","run"))
## horrible hack (but doesn't help); subsample
## w <- which(mL$model=="random")
## mLw <- mL[-sample(w,size=length(w)*9/10,replace=FALSE),]
w <- with(mL,which(model=="random" & variable=="eq_vir"))
rval <- mean(mL$value[w])
mLw <- droplevels(mL[-w,])
mLw$variable <- fixfac2(mLw$variable)

gg_univ <- ggplot(mLw,aes(value,model,fill=model))+
    geom_violinh(width=1)+
    facet_wrap(~variable,scale="free_x")+
    guides(fill=guide_legend(reverse=TRUE))+
    labs(y="",x="")+
    geom_point(data=data.frame(model="random",
                               variable="equilibrium virulence",
                               value=rval),
               pch=22,size=8,show.legend=FALSE) +
    theme(legend.position="none")+
    zero_x_margin

tikz(file="univ.tikz",w
     ## width=12,height=10,
     ## width=10,height=8,
     width=9,height=9*5/6,
     standAlone=TRUE)
print(gg_univ)
dev.off()
