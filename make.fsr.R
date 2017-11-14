#!/usr/bin/env Rscript
## read the data and estimate FSR using ASHR

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

trait.factors <- read_tsv('result/ukbb-fqtl-traits.txt.gz')

## LFSR estimation
trait.factors <- trait.factors %>%
    mutate(pip = 1/(1+exp(-lodds))) %>%
    mutate(pos.prob = pnorm(0, mean = theta, sd = theta.se)) %>%
    mutate(neg.prob = 1 - pos.prob) %>%
    mutate(lfsr = 1 - pip * pmax(pos.prob, neg.prob))

## Calibrate overall FSR rate
trait.factors %>% filter(pip > 0.5) %>%
    summarize(lfsr = mean(lfsr))

