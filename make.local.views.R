#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

## generate local views
in.dir <- argv[1]
out.dir <- argv[2]

options(stringsAsFactors = FALSE)
source('util.R')
library(readr)
library(dplyr)
library(tidyr)
library(scales)

snp.file <- in.dir %&&% '/fqtl.snp-factor.gz'
zscore.file <- in.dir %&&% '/fqtl.zscore.gz'
trait.file <- in.dir %&&% '/fqtl.trait-factor.gz'

## visualize strong SNPs only 
snp.tab <- read_tsv(snp.file)
if(nrow(snp.tab) < 1) q()
zscore.tab <- read_tsv(zscore.file)
if(nrow(zscore.tab) < 1) q()
trait.tab <- read_tsv(trait.file)
if(nrow(trait.tab) < 1) q()

lodds.cutoff <- 2.18 # 0.9

################################################################
trait.cluster.file <- 'result/ukbb-fqtl-traits-slim.txt.gz'

trait.cluster <- read_tsv(trait.cluster.file) %>%
    group_by(CHR, LB, UB, factor, k.sorted) %>%
    slice(which.max(lodds)) %>%
    as.data.frame()

gwas.melt <- zscore.tab %>%
    gather(key = 'trait', value = 'zscore', one_of(unique(trait.tab$trait))) %>%
    mutate(pval = 2 * pnorm(abs(zscore), lower.tail = FALSE)) %>%
    mutate(ln.p = pmin(-log10(pval), 30))

.gwas.kb <- function() {
    function(x) format(x/1e3, big.mark=',')
}

gwas.range <- range(gwas.melt$POS)

.gwas.x.scale <- scale_x_continuous(limits = gwas.range + c(-1000, 1000), expand = c(0, 0),
                                    labels = .gwas.kb(), position = 'bottom')

.gwas.x.scale.flip <- scale_x_continuous(limits = gwas.range + c(-1000, 1000), expand = c(0, 0),
                                         labels = .gwas.kb(), position = 'top')

factors <-
    snp.tab %>% group_by(CHR, LB, UB, factor) %>%
    summarize(snp.lb = min(POS), snp.ub = max(POS)) %>%
    left_join(trait.cluster) %>% na.omit()

if(nrow(factors) < 1) {
    log.msg('There is no significant factor\n\n')
    q()
}

active.traits <- factors %>% select(-trait) %>%
    right_join(trait.tab) %>%
    filter(lodds > 0)

viz.snps <- snp.tab %>% filter(lodds > lodds.cutoff) %>% select(SNP) %>% unique()

################################################################
## show factor by factor
for(k in 1:nrow(factors)) {

    chr <- factors[k, ]$CHR
    ld.idx <- paste(factors[k, 1:3], collapse = '_')
    f <- factors[k, ]$factor
    k.glob <- factors[k, ]$k.sorted
    best.trait <- factors[k, ]$trait

    out.file <- sprintf('%s/K%d/K%d_chr%s_%s.pdf', out.dir, k.glob, k.glob, ld.idx, best.trait)
    dir.create(dirname(out.file), recursive = TRUE)

    if(!file.exists(out.file)) {

        active.traits.sub <- active.traits %>% filter(factor == f)

        ## selected traits
        gwas.active <- gwas.melt %>% filter(trait %in% active.traits.sub$trait)
        snps.active <- snp.tab %>% left_join(gwas.active) %>% na.omit()

        gwas.inactive <- gwas.melt %>% filter(!(trait %in% active.traits.sub$trait))
        snps.inactive <- snp.tab %>% left_join(gwas.inactive) %>% na.omit()

        snps.df <- rbind(snps.active %>% mutate(state = 'active'),
                         snps.inactive %>% mutate(state = 'inactive')) %>%
                             filter(SNP %in% viz.snps$SNP)

        ## strong -> weak traits
        trait.order <- trait.tab %>% filter(factor == f) %>% arrange(lodds) %>% select(trait) %>% as.data.frame()
        snps.df$trait <- factor(snps.df$trait, trait.order$trait)
        print(trait.order)

        gwas.active$trait <- factor(gwas.active$trait, rev(trait.order$trait))
        snps.active$trait <- factor(snps.active$trait, rev(trait.order$trait))
        gwas.inactive$trait <- factor(gwas.inactive$trait, rev(trait.order$trait))

        ## map genomic location to heatmap location
        snps <- snps.df %>% select(SNP, POS) %>% unique() %>% arrange(POS)
        n.snps <- nrow(snps)
        gwas.len <- gwas.range[2] - gwas.range[1]
        gwas.inter <- gwas.len/(n.snps)
        rs.loc <- seq(gwas.range[1] + gwas.inter/2, gwas.range[2] - gwas.inter/2 + 1, by = gwas.inter)
        snps <- cbind(snps, rs.loc)

        snps.df$SNP <- factor(snps.df$SNP, snps$SNP)
        snps.df$trait <- factor(snps.df$trait, trait.order$trait)

        ## show GWAS in inactive traits
        gwas.inactive.med <- gwas.inactive %>%
            group_by(POS) %>%
                summarize(ln.p = median(ln.p))

        ## GWAS manhattan
        p1 <- gg.plot(gwas.active, aes(x = POS, y = ln.p)) +
            geom_vline(data = snps, aes(xintercept = POS), color = 'green', lty = 'solid', size = .5) +
                .gwas.x.scale.flip

        p1 <- p1 +
            geom_point(data = snps.active %>% filter(SNP %in% snps$SNP), aes(size = lodds),
                       pch = 21, fill = 'white', color = 'red') +
                           geom_line(data = gwas.inactive.med, aes(x = POS, y = ln.p), color = 'orange')

        p1 <- p1 +
            geom_point(size = .8, color = 'gray20') +
                facet_wrap(~trait, ncol = 1, scales='free', dir = 'v', strip.position = 'bottom') +
                    scale_size_continuous(range = c(0, 2)) +
                        xlab('genomic location (kb)') + ylab('GWAS (-log10 P)') +
                            theme(legend.position = 'none',
                                  axis.title.x = element_blank(), axis.text = element_text(size = 6))

        ## bridge
        p.12 <-
            gg.plot(snps, aes(x = POS, xend = rs.loc, y = 1, yend = 0)) +
                theme_void() +
                    geom_segment(color = 'gray') +
                        .gwas.x.scale

        ## heatmap of z-scores
        heat.aes <- aes(x = SNP, y = trait, fill = pmin(pmax(-6, zscore), 6))
        scale.fill <- scale_fill_gradient2('z-score', low = '#0000FF', high = '#FF0000')


        p.2 <-
            gg.plot(snps.df %>% filter(state == 'active'), heat.aes) +
                geom_tile(color = 'gray') +
                    geom_text(data = snps.df %>% filter(zscore > 3, state == 'active'),
                              aes(label = signif(zscore, 1)),
                              size = 3, color = 'black')

        p.2 <- p.2 +
            geom_text(data = snps.df %>% filter(zscore < -3, state == 'active'),
                      aes(label = signif(zscore, 1)),
                      size = 3, color = 'white')

        p.2 <- p.2 +
            scale_x_discrete(position = 'top') +
                scale.fill +
                    theme(axis.text.x.top = element_text(angle = 90, vjust = .5, hjust = 0, size = 8),
                          legend.position = 'none', axis.title = element_blank(),
                          axis.text.y = element_text(size=6))

        n.active.trait <- snps.active %>% select(trait) %>% unique() %>% nrow()
        .heights <- c(1 + n.active.trait * 1, .5, 1 + n.active.trait * .25)

        gg <- grid.vcat(list(p1, p.12, p.2), heights = .heights)

        w <- 2 + nrow(snps) * .2
        h <- sum(.heights)

        ggsave(filename = out.file, plot = gg, width = w, height = h, units = 'in', limitsize = FALSE)
    }
}
