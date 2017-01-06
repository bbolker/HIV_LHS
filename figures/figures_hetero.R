library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)
load("../simdata/ev_hetero.rda")
load("../simdata/ev_bd.rda")

zero_x_margin <-
	theme(panel.spacing.x=grid::unit(0, "lines"))

tvec <- c(1:1000)
div <- seq(1, 7, 2)
model <- c("full+epc", "full", "instswitch+epc", "instswitch")

sort_matrix <- function(mat, name){
	list(hetero = data.frame(tvec, mat[,div]), 
			 homo = data.frame(tvec, mat[,-div])) %>%
		lapply(function(x) setNames(x, c("tvec", model))) %>%
		setNames(name) %>%
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

## homo vs hetero

g_base %+% (vir_matFull %>% sort_matrix(c("hetero", "homo")))

g_base %+% (I_matFull %>% sort_matrix(c("hetero", "homo")))

## natural birth and death

g_base %+% (vir_mat %>% sort_matrix(c("no_bd", "bd")))

g_base %+% (I_mat %>% sort_matrix(c("no_bd", "bd")))
