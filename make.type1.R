#!/usr/bin/env Rscript

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

## compare null distributions with observed distributions
null.tab <- read_tsv('result/ukbb-null-kron.txt.gz', col_types = 'iiiccdddd') %>%
    mutate(pip = 1/(1+exp(-lodds)))
                     
n.test <- nrow(null.tab %>% filter(factor != 'total'))

## pip cutoff vs frequency under the null
take.pip.cutoff <- function(.pip) {
    ret <- null.tab %>% filter(factor != 'total') %>%
        summarize(pip.cutoff = .pip, freq = (sum(pip >= .pip) + 1)/(n.test + 1))
    return(ret %>% as.data.frame())
}

pip.cutoff <- 1/(1+exp(-c(seq(-8, 0, by = 0.5))))

pip.type1.tab <- do.call(rbind, lapply(pip.cutoff, take.pip.cutoff)) %>%
    mutate(null = 'correlated')

## another null model
null.tab <- read_tsv('result/ukbb-null.txt.gz', col_types = 'iiiccdddd') %>%
    mutate(pip = 1/(1+exp(-lodds)))

pip.type1.tab.indep <- do.call(rbind, lapply(pip.cutoff, take.pip.cutoff)) %>%
    mutate(null = 'independent')

ret <- rbind(pip.type1.tab, pip.type1.tab.indep)

library(ggplot2)
library(scales)

logit.trans <- trans_new(".logit", function(x) log(x) - log(1-x), function(y) 1/(1+exp(-y)))
l10.trans <- trans_new(".ln10", function(x) pmax(log10(x), -8), function(y) 10^(y))

plt <- gg.plot(ret, aes(x = pip.cutoff, y = freq, color = null, shape = null)) +
    geom_point() + geom_line() +
    scale_x_continuous(breaks = signif(1/(1 + exp(-seq(-6, 0))), 2), trans = logit.trans) +
    scale_y_continuous(breaks = 10^(seq(-8, 0)), trans = l10.trans) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
          legend.title = element_blank(),
          legend.position = c(.99,.99), legend.justification = c(1,1)) +
    ylab('Estimated type I error') + xlab('posterior inclusion probability cutoff')
    
ggsave(filename = 'result/fig_type1_pip.pdf', plot = plt, width = 4, height = 4, units = 'in')

write_tsv(ret, path = gzfile('result/ukbb-null-pip.txt.gz'))
