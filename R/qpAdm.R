#' Perform a 'pre-screening' of the 'right' populations before
#' a qpAdm analysis based on Harney et al., 2020 (bioRxiv)
#'
#' @param data EIGENSTRAT dataset
#' @param candidates Character vector with potential 'outgroup' populations
#' @param left Character vector with target and source populations
#'
#' @return Data frame of all combinations of f4(L_i, L_j; R_k, R_l),
#'   ordered by Zscore of the f4 statistic
#'
#' @export
qpAdm_prescreen <- function(data, candidates, left, Zcutoff = 2) {
  leftcomb <- t(combn(left, 2))
  rightcomb <- t(combn(candidates, 2))

  i <- 1
  quartets <- list()
  for (l in 1:nrow(leftcomb)) {
    for (r in 1:nrow(rightcomb)) {
      quartets[[i]] <- c(leftcomb[l, ], rightcomb[r, ])
      i <- i + 1
    }
  }

  result <- f4(data, quartets = quartets) %>% dplyr::arrange(abs(Zscore))

  right <- dplyr::filter(result, abs(Zscore) > Zcutoff) %>%
    .[, c("Y", "Z")] %>%
    as.matrix %>%
    as.vector %>%
    unique

  list(outgroups = right, screening = result)
}


# Check that the provided object is of the required type
check_type <- function(x, type) {
    if (!inherits(x, type)) {
        stop(glue::glue("Object is not of the type {type}"), call. = FALSE)
    }
}


#' Fit qpAdm models based on the rotation strategy described in
#' Harney et al. 2020 (bioRxiv)
#'
#' @param data EIGENSTRAT dataset
#' @param target Target population that is modeled as admixed
#' @param candidates Potential candidates for sources and outgroups
#' @param nsources Number of sources to pull from the candidates
#' @param ncores Number of CPU cores to utilize for model fitting
#'
#' @return qpAdm list with proportions, ranks and subsets elements (as
#'     with a traditional qpAdm run)
#'
#' @export
qpAdm_rotation <- function(data, target, candidates, minimize = FALSE, nsources = 2, ncores = 1) {
    check_type(data, "EIGENSTRAT")

    ## generate combinations of possible sources and outgroups
    sources <- t(combn(candidates, nsources))
    sources_outgroups <- unlist(lapply(1:nrow(sources), function(i) {
        outgroups <- setdiff(candidates, sources[i, ])
        if (minimize) {
            outgroups <- unlist(lapply((nsources + 1):length(outgroups), function(nout) {
                outcomb <- t(combn(outgroups, nout))
                lapply(1:nrow(outcomb), function(j) outcomb[j, ])
            }), recursive = FALSE)
        } else {
            outgroups <- list(outgroups)
        }

        lapply(outgroups, function(out) { list(sources = sources[i, ], outgroups = out) })
    }), recursive = FALSE)

    ## run qpAdm for all combinations of sources and outgroups
    results_list <- parallel::mclapply(sources_outgroups, function(x) {
        result <- qpAdm(
            data,
            target = target, sources = x$sources, outgroups = x$outgroups
        )

        names(x$sources) <- paste0("source", 1:nsources)
        sources_df <- as.data.frame(t(as.matrix(x$sources)))

        ## rename sources and stderr columns (by default, proportion and
        ## stderr columns are named based on the source populations - we
        ## don't want that here because we want to merge all the individual
        ## proportion tables, columns have to have the same name)
        names(result$proportions)[2:(1 + nsources)] <- paste0("prop", 1:nsources)
        names(result$proportions)[4:(3 + nsources)] <- paste0("stderr", 1:nsources)
        ## add source names as two new columns
        result$proportions <- cbind(result$proportions, sources_df) %>%
            dplyr::mutate(outgroups = paste0(x$outgroups, collapse = " & "))
        ## rearrange columns
        result$proportions <- dplyr::select(
            result$proportions, target, names(x$sources), outgroups, pvalue,
            dplyr::everything()
        )

        ## reformat rank table
        names(result$subsets)[7:(6 + nsources)] <- paste0("prop", 1:nsources)
        ## add source names as two new columns
        result$subsets <- cbind(result$subsets, sources_df)
        result$subsets <- dplyr::select(
            result$subsets, target, names(x$sources), pattern,
            dplyr::everything()
        )

        result
    }, mc.cores = ncores)
    
    ## extract log information before further processing of the results
    log_lines <- sapply(results_list, function(i) attr(i, "log_output"))
    names(log_lines) <- paste0("m", seq_along(sources_outgroups))

    proportions <- dplyr::bind_rows(lapply(results_list, `[[`, "proportions"))
    ranks <- dplyr::bind_rows(lapply(results_list, `[[`, "ranks"))
    subsets <- dplyr::bind_rows(lapply(results_list, `[[`, "subsets"))

    ## add model identifier to each row in the proportions table, ...
    models <- paste0("m", seq_along(sources_outgroups))
    proportions$model <- models
    proportions <- dplyr::as_tibble(proportions) %>% dplyr::select(model, everything())
    ## ... ranks table, ...
    ranks$model <- sort(rep(models, 2))
    ranks <- dplyr::as_tibble(ranks) %>% dplyr::select(model, everything())
    ## and subsets table
    subsets$model <- sort(rep(models, 1 + 2^(nsources - 1)))
    subsets <- dplyr::as_tibble(subsets) %>% dplyr::select(model, everything())

    ## bind all tables together and add log information
    results <- list(proportions = proportions, ranks = ranks, subsets = subsets)
    attr(results, "command") <- "qpAdm_rotation"
    attr(results, "log_output") <- log_lines
    class(results) <- c("admixr_result", class(results))

    results
}


#' Filter qpAdm rotation results for only 'sensible' models
#'
#' Filter for p-value larger than a specified cuttof and admixture
#' proportions between 0 and 1.
#'
#' @param x Output of a qpAdm_rotation() function
#' @param p p-value cutoff (default 0: will only filter for sensible
#'     admixture proportions)
#'
#' @return qpAdm_rotation object filtered down based on p-value
qpAdm_filter <- function(x, p = 0) {
    ## get positions of columns with estimated admixture proportions
    prop_columns <- stringr::str_which(names(x$proportions), "prop")
    
    ## find out rows/models for which all proportions are in [0, 1] and
    ## pvalue is larger than the required cutoff
    pvalue <- x$proportions$pvalue > p
    constr_0 <- x$proportions[, prop_columns] >= 0
    constr_1 <- x$proportions[, prop_columns] <= 1
    constr <- apply(pvalue & constr_0 & constr_1, 1, all)

    ## filter all three sub-tables to only those models that fit the criteria
    x$proportions <- dplyr::arrange(x$proportions[constr, ], -pvalue)
    x$ranks <- x$ranks[x$ranks$model %in% x$proportions$model, ]
    x$subsets <- x$subsets[x$subsets$model %in% x$proportions$model, ]

    ## filter also only to relevant remaining log output information
    attr(x, "log_output") <- attr(x, "log_output")[unique(x$proportions$model)]

    x
}
