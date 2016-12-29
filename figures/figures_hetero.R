library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)
load("../simdata/ev_hetero.rda")

zero_x_margin <-
	theme(panel.spacing.x=grid::unit(0, "lines"))

tvec <- c(1:1000)
het <- seq(1, 7, 2)
model <- c("full+epc", "full", "instswitch+epc", "instswitch")

sort_matrix <- function(mat){
	list(hetero = data.frame(tvec, mat[,het]), 
			 homo = data.frame(tvec, mat[,-het])) %>%
		lapply(function(x) setNames(x, c("tvec", model))) %>%
		setNames(c("hetero", "homo")) %>%
		bind_rows(.id = "type") %>%
		as_data_frame %>%
		gather(key, value, -type, -tvec) %>%
		mutate(key = factor(.$key, levels = model)) %>%
		filter(tvec < 600)
}

g_base <- ggplot(data = NULL, aes(tvec, value, col = type)) +
	geom_line() +
	scale_x_continuous(expand = c(0.005,0), breaks = seq(0, 500, 100), labels = c("0", "", "200", "", "400", "")) +
	facet_grid(.~key) +
	zero_x_margin +
	theme(panel.grid = element_blank())

g_base %+% (vir_matFull %>% sort_matrix)

g_base %+% (I_matFull %>% sort_matrix)
