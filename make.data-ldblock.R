#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

stat.file <- argv[1] # e.g., 'UKBB/blood_EOSINOPHIL_COUNT.sumstats.gz'
ld.file <- argv[2] # e.g., 'ldblocks/EUR/fourier_ls-all.bed'
out.dir <- argv[3]

options(stringsAsFactors = FALSE)
library(readr)
library(dplyr)
source('util.R')

trait <- basename(stat.file) %>%
  gsub(pattern = '.sumstats.gz', replacement = '')

ld.tab <- read_tsv(ld.file) %>%
    rename(CHR = chr, START = start, END = stop) %>%
    mutate(CHR = gsub(CHR, pattern = 'chr', replacement = ''))

stat.tab <- read_tsv(stat.file)

for(rr in 1:nrow(ld.tab)) {

    chr <- ld.tab$CHR[rr]
    ss <- ld.tab$START[rr]
    ee <- ld.tab$END[rr]

    out.file <- out.dir %&&% '/' %&&% chr %&&% '/' %&&% rr %&&% '/' %&&% trait %&&% '.stat.gz'
    dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

    if(!file.exists(out.file)) {
        .temp <- stat.tab %>%
            filter(CHR == chr, POS >= ss, POS < ee, EAF >= .05, INFO >= .9) %>%
                mutate(EAF = signif(EAF,2), P = signif(P, 2))
        
        write_tsv(.temp, path = gzfile(out.file))
        rm(.temp)
        gc()
        log.msg('Wrote %s file\n', out.file)
    }
}

log.msg('Successfully broke down data\n\n')
