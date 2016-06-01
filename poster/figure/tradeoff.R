source("fig_setup.R")

pars2 <- HIVpars.base
pars2$alphaDist <- c(min=0,max=9,delta=0.25)
pars2 <- expand(pars2)
dd <-with(pars2,
     melt(data.frame(alpha=alpha,"transmission probability"=Beta2,
                     "duration (years)"=Duration2,check.names=FALSE),
          id.vars="alpha"))
gg_tradeoff <- ggplot(dd,aes(alpha,value))+facet_wrap(~variable,scale="free")+
    geom_line(size=7)+
    labs(x= Lspvl,y="")

tikz(file="tradeoff.tikz",width=12,height=8,standAlone=TRUE)
print(gg_tradeoff)
dev.off()
