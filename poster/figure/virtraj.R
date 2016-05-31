source("fig_setup.R")

ss <- subset(sum_list[["vir_mat"]],tvec<750)
ss$model <- factor(ss$model, c("random",
         "pairform+epc", "pairform","instswitch+epc","instswitch","implicit"))

gg_virtraj <- ggplot(ss,
                     aes(tvec,mean,ymin=lwr,ymax=upr))+
    geom_line(aes(colour=model,linetype=model),lwd=1)+
    geom_ribbon(aes(fill = model), alpha=0.3)+
    labs(x="time (years)",y="population mean set-point viral load ($\\log_{10}$)") +
    theme(axis.title.y = element_text(margin = margin(0,10,0,0)),
          axis.title.x = element_text(margin = margin(10,0,0,0))) +
    facet_wrap(~model, nrow = 3) +
    zero_margin+
    theme(legend.position="none")+
    scale_x_continuous(breaks=c(0,600))

tikz(file="virtraj.tikz",width=8,height=11,standAlone=TRUE)
print(gg_virtraj)
dev.off()
