#!/usr/bin/env Rscript
## Identify clusters of tissue patterns

library(readr)
library(dplyr)
library(tidyr)
source('util.R')

library(ClusterR)

trait.tab <- read_tsv('result/ukbb-fqtl-traits.txt.gz')
snp.tab <- read_tsv('result/ukbb-fqtl-snps.txt.gz')
trait.out.file <- 'result/ukbb-fqtl-traits-slim.txt.gz'
snp.out.file <- 'result/ukbb-fqtl-snps-slim.txt.gz'
centroid.out.file <- 'result/ukbb-centroid-slim.txt.gz'
pair.out.file <- 'result/ukbb-trait-pairs.txt.gz'

lodds.cutoff <- 2
## only include factors, which contains at least one lodds >= 2
data.idx <- 
    trait.tab %>% filter(lodds >= lodds.cutoff) %>%
    select(CHR, LB, UB, trait, factor) %>%
    unique()

data.tab <- trait.tab %>% full_join(data.idx) %>%
    na.omit() %>% 
    select(CHR, LB, UB, trait, factor, lodds) %>%
    mutate(lodds = pmin(pmax(lodds, -3), 3)) %>%
    spread(key = trait, value = lodds)

M <- data.tab[, -(1:4)] %>% as.matrix()

################################################################
## run clustering

opt <- Optimal_Clusters_KMeans(M, max_clusters = 100,
                               num_init = 7,
                               plot_clusters = FALSE,
                               verbose = TRUE,
                               criterion = 'BIC',
                               seed = 13)

K <- which.min(opt)

clust.out <- KMeans_rcpp(M, clusters = K, num_init = 50, max_iters = 200,
                         verbose = TRUE, seed = 13)

################################################################
## reorder centroids
C <- clust.out$centroid
ko <- row.order(C > 0)    # cluster order
to <- row.order(t(C) > 0) # trait order
traits <- colnames(M)
traits.ordered <- traits[to]

## pair statistics
.take.t1.t2 <- function(t1, t2, .tab = data.tab, lodds.cutoff = 0) {
    .temp <- .tab %>% select_(.dots = c(t1, t2))
    n1 <- sum(.temp[, 1] > lodds.cutoff)
    n2 <- sum(.temp[, 2] > lodds.cutoff)
    n12 <- sum(.temp[, 2] > lodds.cutoff & .temp[, 1] > lodds.cutoff)

    ret <- data.frame(t1, t2, jaccard = n12 / pmax(n1 + n2 - n12, 1),
                      p1.2 = n12 / pmax(n1, 1),
                      p2.1 = n12 / pmax(n2, 1))
    return(ret)
}

.func.j <- function(j) {
    .take.t1.t2(trait.pairs[1, j], trait.pairs[2, j], .tab = data.tab, lodds.cutoff = 0)
}

trait.pairs <- combn(traits.ordered, 2)

pair.stat.tab <- do.call(rbind, lapply(1:ncol(trait.pairs), .func.j))

pair.stat <- rbind(data.frame(t1 = pair.stat.tab$t1, t2 = pair.stat.tab$t2, p.yx=pair.stat.tab$p2.1),
                   data.frame(t1 = pair.stat.tab$t2, t2 = pair.stat.tab$t1, p.yx=pair.stat.tab$p1.2))

write_tsv(pair.stat, path = gzfile(pair.out.file))

################################################################
## output trait clusters with membership assigned

clust.idx <- data.frame(k = clust.out$clusters) %>%
    mutate(data.pos = 1:n()) %>%
    left_join(data.frame(k = ko) %>% mutate(k.sorted = 1:n())) %>%
    select(data.pos, k.sorted)

data.tab.cluster <- data.tab %>% mutate(data.pos = 1:n()) %>% left_join(clust.idx) %>%
    select(-data.pos)

trait.tab.slim <- data.tab.cluster %>%
    gather(key = trait, value = lodds, traits) %>%
    filter(lodds > 0)

write_tsv(trait.tab.slim, path = gzfile(trait.out.file))

################################################################
## output snp clusters with membeship assignment

snp.tab.slim <- snp.tab %>%
    left_join(trait.tab.slim %>% select(CHR, LB, UB, factor, k.sorted) %>% unique()) %>%
              na.omit()

write_tsv(snp.tab.slim, path = gzfile(snp.out.file))

## cluster-specific stats
.snps <- snp.tab.slim %>% select(SNP, k.sorted) %>% unique() %>%
    group_by(k.sorted) %>% summarize(num.snps = n())

.factors <- trait.tab.slim %>% select(CHR, LB, UB, factor, k.sorted) %>%
    mutate(sz = UB - LB) %>%
    unique() %>% group_by(k.sorted) %>%
    summarize(num.factors = n(), size = sum(sz))

clust.stat <- left_join(.snps, .factors) %>% rename(cluster = k.sorted)


## trait-specific stats
.snps <- snp.tab.slim %>% select(CHR, LB, UB, factor, SNP) %>%
    left_join(trait.tab.slim %>% select(CHR, LB, UB, factor, trait)) %>%
    select(trait, SNP) %>% unique() %>%
    group_by(trait) %>%
    summarize(num.snps = n())

.factors <- snp.tab.slim %>% select(CHR, LB, UB, factor, SNP) %>%
    left_join(trait.tab.slim %>% select(CHR, LB, UB, factor, trait)) %>%
    select(CHR, LB, UB, trait) %>% unique() %>%
    group_by(trait) %>%
    summarize(num.factors = n())

trait.stat <- left_join(.snps, .factors)

################################################################
## Draw centroids

colnames(clust.out$centroid) <- traits

centroid.tab <- cbind(cluster = 1:K, clust.out$centroid) %>% as.data.frame() %>%
    gather(key = trait, value = lodds, traits) %>%
    mutate(pip = 1/(1+exp(-lodds)))

centroid.tab$trait <- factor(centroid.tab$trait, traits.ordered)
centroid.tab$cluster <- factor(centroid.tab$cluster, ko, 1:K)

write_tsv(centroid.tab, path = gzfile(centroid.out.file))

plt <-
    gg.plot(centroid.tab, aes(x = cluster, y = trait, fill = pip)) +
    geom_tile(color = 'gray') + xlab('pleiotropic genomic regions') + ylab('UKBB 45 traits') +
    scale_fill_gradientn('PIP', colors = c('#FFFFFF', '#FF0000'),
                         breaks = c(0, .05, .25, .5, .75, .95, 1))

plt <- plt +
    theme(legend.position = c(.98, .98), legend.justification = c(1,1),
          legend.background = element_rect(fill = 'white', color = 'white'),
          axis.text.x = element_text(size = 4))

trait.stat$trait <- factor(trait.stat$trait, traits.ordered)

p.0 <- gg.plot() + geom_blank() + theme_void()

clust.stat$cluster <- factor(clust.stat$cluster, 1:K)

p.1a <- 
    gg.plot(clust.stat)+
    geom_segment(aes(x = cluster, xend = cluster, y = 0, yend = num.factors), size = 2, color = 'gray')+
    geom_text(aes(x = cluster, y = num.factors + 0.1, label = num.factors), angle = 90, vjust = .5, hjust = 0, size = 3)+
    ylab('# LD blocks') + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

p.1b <- 
    gg.plot(clust.stat)+
    geom_segment(aes(x = cluster, xend = cluster, y = 0, yend = size), size = 2, color = 'gray')+
    ylab('# total SNPs') + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

p.1c <-
    gg.plot(clust.stat)+
    geom_segment(aes(x = cluster, xend = cluster, y = 0, yend = num.snps), size = 2, color = 'gray')+
    geom_text(aes(x = cluster, y = num.snps + 10, label = num.snps), angle = 90, vjust = .5, hjust = 0, size = 3)+
    ylab('# SNPs') + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

p.2 <-
    gg.plot(trait.stat, aes(x = 0, xend = num.snps, y = trait, yend = trait)) +
    geom_segment(color = 'gray', size = 2) +
    geom_text(aes(x = pmax(num.snps, 500), label = num.snps), size = 3, hjust = 1) +
    xlab('# SNPs') + theme(axis.text.y = element_blank(), axis.title.y = element_blank())

g.temp <- match.widths(list(p.1a, p.1b, p.1c, plt))

## g.top <- match.heights.grob(c(g.temp[1], list(ggplotGrob(p.0))))

g.bottom <- match.heights.grob(c(g.temp[4], list(ggplotGrob(p.2))))

g.top <- c(g.temp[1], list(ggplotGrob(p.0)),
           g.temp[2], list(ggplotGrob(p.0)),
           g.temp[3], list(ggplotGrob(p.0)))

gg <-
    grid.arrange(grobs = c(g.top, g.bottom), ncol = 2, newpage = TRUE,
             widths = c(9, 1), heights = c(2, 2, 2, 10))

ggsave(filename = 'result/fig_trait_clusters.pdf',
       plot = gg,
       width = K * .1 + 5,
       height = 6 + length(traits.ordered) * .1,
       limitsize = FALSE,
       units = 'in')


################################################################
## Draw trait-trait sharing -- bidirectional posterior probability matrix

pair.df <- rbind(data.frame(t1 = traits, t2 = traits, p.yx = 1),
                 pair.stat)
pair.df$t1 <- factor(pair.df$t1, rev(traits.ordered))
pair.df$t2 <- factor(pair.df$t2, traits.ordered)

plt2 <-
    gg.plot(pair.df, aes(x = t1, y = t2, fill = p.yx)) +
    geom_tile(color = 'black') +
    geom_text(data = pair.df %>% filter(t1 != t2, p.yx >= .2), aes(label = floor(10*p.yx)),
              size = 3) + 
    scale_x_discrete(position = 'top') +
    scale_fill_gradientn('P(row|col)', colors = c('white', 'white', 'orange', 'red'),
                         breaks = c(0.1, 0.25, 0.5, 0.8, 1),
                         trans = 'sqrt') +
    theme(axis.text.x.top = element_text(angle=60, vjust = 0, hjust=0),
          axis.title = element_blank())
          
ggsave(filename = 'result/fig_trait_conditional.pdf',
       plot = plt2,
       width = 4 + length(traits.ordered) * .15,
       height = 3 + length(traits.ordered) * .15,
       limitsize = FALSE,
       units = 'in')

log.msg('Done\n')
