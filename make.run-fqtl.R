#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

in.dir <- argv[1]               # e.g., in.dir = '/broad/hptmp/ypp/ukbb/tempdata/1/90/'
plink.hdr <- argv[2]            # e.g., plink.hdr = '1KG/plink/chr1'
N.TRAITS <- as.integer(argv[3]) # e.g., N.TRAITS = 45
out.hdr <- argv[4]              # e.g., out.hdr = 'output'

options(stringsAsFactors = FALSE)
source('util.R')

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
snp.factor.file <- out.hdr %&&% '.snp-factor.gz'
trait.factor.file <- out.hdr %&&% '.trait-factor.gz'
zscore.file <- out.hdr %&&% '.zscore.gz'
var.file <- out.hdr %&&% '.var.gz'

.files <- c(snp.factor.file, trait.factor.file, zscore.file, var.file)

if(all(sapply(.files, file.exists))) {
    log.msg('Files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

library(readr)
library(dplyr)
library(tidyr)
library(zqtl)

## Take European ancestry
ref.sample.info <- read_tsv('20130606_g1k.ped') %>%
    filter(Population %in% c('CEU', 'TSI', 'FIN', 'GBR', 'IBS'))

stat.files <- list.files(path = in.dir, pattern = '.stat.gz', full.names = TRUE)

if(length(stat.files) < N.TRAITS) {
    log.msg('Fewer number of stat files: %d < %d\n\n',
            length(stat.files), N.TRAITS)
    q()
}

read.stat <- function(stat.file) {
    trait <- basename(stat.file) %>%
        gsub(pattern='.stat.gz', replacement = '')

    ret <- read_tsv(stat.file) %>%
        group_by(CHR, SNP, POS, A1, A2, REF) %>%
            summarize(Beta = mean(Beta), se = mean(se)) %>%
                mutate(trait = trait)
    return(ret)
}

## plink data
plink <- read.plink(plink.hdr)
colnames(plink$BIM) <- c('CHR', 'SNP', '.', 'POS', 'plink.A1', 'plink.A2')
plink.snps <- plink$BIM %>% mutate(plink.pos = 1:n()) %>%
    select(CHR, SNP, POS, plink.A1, plink.A2, plink.pos)

## match directionality between 1kg and ukbb
stat.tab <- do.call(rbind, lapply(stat.files, read.stat)) %>%
    left_join(plink.snps) %>% na.omit() %>%
        filter((A1 == plink.A2 & A2 == plink.A1) | (A1 == plink.A1 & A2 == plink.A2)) %>%
            mutate(Beta = ifelse(A1 == plink.A1, -Beta, Beta)) %>%
                as.data.frame()

nsnps.traits <- stat.tab %>% select(trait, POS) %>% unique() %>%
    group_by(trait) %>% summarize(n = n())

## check if there were different number of SNPs
max.nsnps <- max(nsnps.traits$n)
discrepancy <- nsnps.traits %>% filter(n < (max.nsnps/2))
traits <- stat.tab$trait %>% unique() %>% sort()

if(discrepancy %>% nrow() > 0 || length(traits) < N.TRAITS) {
    library(pander)
    log.msg('Some traits might have preprocessing errors:\n%s\n%s\n\n',
            pandoc.table.return(nsnps.traits, style = 'simple'),
            paste(stat.files, collapse = '\n'))
    q()
}

## take beta and se mat
chr <- stat.tab$CHR[1]

beta.tab <- stat.tab %>% select(CHR, SNP, POS, A1, A2, REF, plink.pos, Beta, trait) %>%
    spread(key = trait, value = Beta) %>% arrange(POS, A1) %>% na.omit() %>%
        as.data.frame()

se.tab <- stat.tab %>% select(CHR, SNP, POS, A1, A2, REF, plink.pos, se, trait) %>%
    spread(key = trait, value = se) %>% arrange(POS, A1) %>% na.omit() %>%
        as.data.frame()

zscore.tab <- stat.tab %>% select(CHR, SNP, POS, A1, A2, REF, plink.pos, Beta, se, trait) %>%
    mutate(z = signif(Beta / se, 2)) %>% select(-Beta, -se) %>%
        spread(key = trait, value = z) %>% arrange(POS, A1) %>% na.omit() %>%
            as.data.frame()

if(nrow(beta.tab) < 2) {
    log.msg('Too few number of SNPs\n\n')
    write_tsv(data.frame(NULL), path = gzfile(var.file))
    write_tsv(data.frame(NULL), path = gzfile(snp.factor.file))
    write_tsv(data.frame(NULL), path = gzfile(trait.factor.file))
    write_tsv(data.frame(NULL), path = gzfile(zscore.file))
    q()
}

rm(stat.tab)
gc()

colnames(plink$FAM) <- c('Individual ID', 'Family ID', '1', '2', '3', '4')

fam <- plink$FAM %>% mutate(fam.pos = 1:n()) %>%
    left_join(ref.sample.info, by = 'Individual ID') %>%
    na.omit()

X <- plink$BED %c% beta.tab$plink.pos %r% fam$fam.pos

plink.snps <- plink.snps %r% beta.tab$plink.pos
rm(plink)
gc()

################################################################
tab2mat <- function(tab) {
    tab %>% select_(.dots = traits) %>% as.matrix()
}

beta.mat <- plink.snps %>% select(CHR, SNP, POS) %>% left_join(beta.tab) %>% na.omit() %>%
    tab2mat()

se.mat <- plink.snps %>% select(CHR, SNP, POS) %>% left_join(se.tab) %>% na.omit() %>%
    tab2mat()

rm(se.tab)
rm(beta.tab)
gc()

zscore.tab <- plink.snps %>% select(CHR, SNP, POS) %>% left_join(zscore.tab) %>% na.omit()
log.msg('Constructed data\n\n')

################################################################
K <- max(min(length(traits), ncol(X)), 1)

vb.opt <- list(tau.lb = -10, tau.ub = -4, pi.lb = -4, pi.ub = -2, do.hyper = TRUE,
               do.stdize = TRUE, eigen.tol = 1e-2, weight.y = FALSE, gammax = 1e4,
               svd.init = TRUE, jitter = 0.01, vbiter = 5000, tol = 1e-8, rate = 1e-2,
               k = K)

z.out <- fit.zqtl(effect = beta.mat, effect.se = se.mat,
                  X = X, factored = TRUE, options = vb.opt)

log.msg('Finished fQTL estimation\n\n')

LD.info <- plink.snps %>% summarize(CHR = min(CHR), LB = min(POS), UB = max(POS))

## Variance of each factor
theta.snp <- z.out$param.left$theta
theta.trait <- z.out$param.right$theta

take.var.k <- function(k) {

    theta.k <- (theta.snp %c% k) %*% t(theta.trait %c% k)
    theta.k <- theta.k * se.mat ## scale by SE
    eta.k <- sweep(z.out$Vt %*% theta.k, 1, sqrt(z.out$D2), `*`)
    var.k <- apply(eta.k^2, 2, sum)

    return(data.frame(trait = as.character(names(var.k)), factor = k, var = var.k))
}

theta.tot <- (theta.snp %*% t(theta.trait)) * se.mat ## scale by SE
eta.tot <- sweep(z.out$Vt %*% theta.tot, 1, sqrt(z.out$D2), `*`)
var.tot <- data.frame(trait = as.character(traits), factor = 'total', var = apply(eta.tot^2, 2, sum)) 

var.tab <- do.call(rbind, lapply(1:K, take.var.k))
rownames(var.tab) <- NULL
var.tab <- rbind(var.tab, var.tot) %>%
    mutate(var = signif(var, 4))

var.tab <- data.frame(LD.info, var.tab)

log.msg('Calculated Variance\n\n')

right.tab <- melt.effect(z.out$param.right, traits, 1:K) %>%
    rename(theta.se = theta.var, trait = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

right.tab <- data.frame(LD.info, right.tab)

left.tab <- melt.effect(z.out$param.left, zscore.tab$SNP, 1:K) %>%
    filter(lodds > -2) %>%
        rename(theta.se = theta.var, SNP = row, factor = col) %>%
            mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
                left_join(zscore.tab %>% select(SNP, POS, A1, A2, REF)) %>%
                    as.data.frame()

if(nrow(left.tab) > 0){
    left.tab <- data.frame(LD.info, left.tab)
} else {
    left.tab <- data.frame(NULL)
}

write_tsv(var.tab, path = gzfile(var.file))
write_tsv(left.tab, path = gzfile(snp.factor.file))
write_tsv(right.tab, path = gzfile(trait.factor.file))
write_tsv(zscore.tab, path = gzfile(zscore.file))

log.msg('Successfully finished fQTL\n\n')
