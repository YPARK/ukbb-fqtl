#!/usr/bin/env Rscript

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

.take.trait.tab <- function(.blk, lodds.cutoff = 0) {

    .trait.file <- .blk %&&% '/fqtl.trait-factor.gz'

    if(!file.exists(.trait.file)) return(NULL)

    if(lodds.cutoff <= -Inf) {
        .trait.tab <- suppressMessages(read_tsv(.trait.file, col_names = TRUE))
        return(.trait.tab %>% as.data.frame())
    }

    .snp.file <- .blk %&&% '/fqtl.snp-factor.gz'

    if(!file.exists(.snp.file)) return(NULL)

    .snp.tab <- suppressMessages(read_tsv(.snp.file, col_names= TRUE) %>%
                                     dplyr::filter(lodds > lodds.cutoff))
    
    if(nrow(.snp.tab) > 0) {        
        .factors <- .snp.tab %>% dplyr::select(factor) %>% unique()        
        .trait.tab <- suppressMessages(read_tsv(.trait.file, col_names = TRUE) %>%
                                           dplyr::filter(factor %in% .factors$factor))
        return(.trait.tab %>% as.data.frame())
    }
    return(NULL)
}

.take.var.tab <- function(.blk) {
    .var.file <- .blk %&&% '/fqtl.var.gz'
    if(!file.exists(.var.file)) return(NULL)
    .var.tab <- suppressMessages(read_tsv(.var.file, col_names = TRUE, col_types = 'iiiccd'))
    return(.var.tab %>% as.data.frame())
}

.take.snp.tab <- function(.blk, lodds.cutoff = 0) {
    .snp.file <- .blk %&&% '/fqtl.snp-factor.gz'
    if(!file.exists(.snp.file)) return(NULL)

    .snp.tab <- suppressMessages(read_tsv(.snp.file, col_names= TRUE) %>%
                                     dplyr::filter(lodds > lodds.cutoff))
    
    if(nrow(.snp.tab) > 0) {        
        return(.snp.tab)
    }
    return(NULL)
}

result.dir <- 'result/fqtl/'

.list.files <- function(chr) list.files(result.dir %&&% chr %&&% '/', full.names = TRUE)

total.blks <- do.call(c, lapply(1:22, .list.files))

trait.tab <- do.call(rbind, lapply(total.blks, .take.trait.tab))
write_tsv(trait.tab, path = gzfile('result/ukbb-fqtl-traits.txt.gz'))
rm(trait.tab); gc();

trait.tot.tab <- do.call(rbind, lapply(total.blks, .take.trait.tab, lodds.cutoff = -Inf))
write_tsv(trait.tot.tab, path = gzfile('result/ukbb-fqtl-tot-traits.txt.gz'))
rm(trait.tot.tab); gc();

var.tot.tab <- do.call(rbind, lapply(total.blks, .take.var.tab))
write_tsv(var.tot.tab, path = gzfile('result/ukbb-fqtl-tot-var.txt.gz'))
rm(var.tot.tab); gc();

snp.tab <- do.call(rbind, lapply(total.blks, .take.snp.tab))
write_tsv(snp.tab, path = gzfile('result/ukbb-fqtl-snps.txt.gz'))

