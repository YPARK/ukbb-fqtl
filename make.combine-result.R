#!/usr/bin/env Rscript

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

.take.trait.tab <- function(.blk, lodds.cutoff = 0) {
    .snp.file <- .blk %&&% '/fqtl.snp-factor.gz'
    if(!file.exists(.snp.file)) return(NULL)

    .snp.tab <- suppressMessages(read_tsv(.snp.file, col_names= TRUE) %>%
                                     dplyr::filter(lodds > lodds.cutoff))
    
    if(nrow(.snp.tab) > 0) {        
        .factors <- .snp.tab %>% dplyr::select(factor) %>% unique()        
        .trait.tab <- suppressMessages(read_tsv(.blk %&&% '/fqtl.trait-factor.gz', col_names = TRUE) %>%
                                           dplyr::filter(factor %in% .factors$factor))
        return(.trait.tab %>% as.data.frame())
    }
    return(NULL)
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

snp.tab <- do.call(rbind, lapply(total.blks, .take.snp.tab))

write_tsv(trait.tab, path = gzfile('result/ukbb-fqtl-traits.txt.gz'))
write_tsv(snp.tab, path = gzfile('result/ukbb-fqtl-snps.txt.gz'))

