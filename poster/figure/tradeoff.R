source("fig_setup.R")

dd <-with(pars,
     melt(data.frame(alpha=alpha,"transmission probability"=Beta2,
                     "duration (years)"=Duration2,check.names=FALSE),
          id.vars="alpha"))
gg_tradeoff <- ggplot(dd,aes(alpha,value))+facet_wrap(~variable,scale="free")+
    geom_line(size=8)+
    labs(x= Lspvl,y="")

tikz(file="tradeoff.tikz",width=12,height=8,standAlone=TRUE)
print(gg_tradeoff)
dev.off()
