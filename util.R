options(stringsAsFactors = FALSE)

log.msg <- function(...) {
    cat('[', date() ,']', sprintf(...), file = stderr())
}

`%&&%` <- function(a, b) paste(a, b, sep = '')

`%c%` <- function(a, b) a[, b, drop = FALSE]

`%r%` <- function(a, b) a[b, , drop = FALSE]

.NA <- function(nrow, ncol) {
    matrix(NA, nrow, ncol)
}

.rnorm <- function(nrow, ncol) {
    matrix(rnorm(nrow * ncol), nrow, ncol)
}

.sample <- function(nrow, ncol, pop = c(1,-1)) {
    matrix(sample(pop, nrow * ncol, TRUE), nrow, ncol)
}

melt.effect <- function(effect.obj, .rnames, .cnames) {

    .melt.effect <- function(mat.list, val.names, i) {
        mat <- signif(mat.list[[i]], 4)
        val <- val.names[[i]]
        colnames(mat) <- .cnames
        ret <- mat %>% as.data.frame() %>% mutate(row = .rnames) %>%
            gather(key = 'col', value = 'value', .cnames) %>%
                rename_(.dots = setNames('value', val)) %>%
                    as.data.frame()
        return(ret)
    }

    melt.list <- lapply(seq_along(effect.obj), .melt.effect,
                        mat.list = effect.obj,
                        val.names = names(effect.obj))

    ret <- Reduce(function(...) left_join(..., by = c('row', 'col')), melt.list) %>%
        as.data.frame()
    return(ret)
}
