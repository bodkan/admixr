#' Fit qpAdm models based on the rotation strategy described in
#' Harney et al. 2020 (bioRxiv)
#'
#' @param data EIGENSTRAT dataset
#' @param target Target population that is modeled as admixed
#' @param candidates Potential candidates for sources and outgroups
#' @param minimize Test also all possible subsets of outgroups? (default TRUE)
#' @param nsources Number of sources to pull from the candidates
#' @param ncores Number of CPU cores to utilize for model fitting
#' @param fulloutput Report also 'ranks' and 'subsets' analysis from
#'     qpAdm in addition to the admixture proportions results? (default FALSE)
#'
#' @return qpAdm list with proportions, ranks and subsets elements (as
#'     with a traditional qpAdm run) or just the proportions
#'     (determined by the value of the 'fulloutput' argument)
#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # find the set of most likely two-source qpAdm models of
#' # a French individual - produce only the 'proportions'
#' # qpAdm summary
#' models <- qpAdm_rotation(
#'     data = snps,
#'     target = "French",
#'     candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija",
#'                    "Altai", "Denisova", "Chimp"),
#'     minimize = TRUE,
#'     nsources = 2,
#'     ncores = 2,
#'     fulloutput = FALSE
#' )
#' }
#'
#' @importFrom utils combn
#' @export
qpAdm_rotation <- function(data, target, candidates, minimize = TRUE, nsources = 2, ncores = 1, fulloutput = FALSE) {
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
            dplyr::mutate(outgroups = paste0(x$outgroups, collapse = " & "),
                          noutgroups = length(x$outgroups))
        ## rearrange columns
        result$proportions <- dplyr::select(
            result$proportions, target, names(x$sources), outgroups, noutgroups, pvalue,
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
    proportions <- dplyr::as_tibble(proportions) %>% dplyr::select(model, dplyr::everything())
    ## ... ranks table, ...
    ranks$model <- sort(rep(models, 2))
    ranks <- dplyr::as_tibble(ranks) %>% dplyr::select(model, dplyr::everything())
    ## and subsets table
    subsets$model <- sort(rep(models, 1 + 2^(nsources - 1)))
    subsets <- dplyr::as_tibble(subsets) %>% dplyr::select(model, dplyr::everything())

    ## add metadata to the results object
    if (fulloutput)
        results <- list(proportions = proportions, ranks = ranks, subsets = subsets)
    else
        results <- proportions
        
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
#'
#' 
#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # find the set of most likely two-source qpAdm models of
#' # a French individual - produce only the 'proportions'
#' # qpAdm summary
#' models <- qpAdm_rotation(
#'     data = snps,
#'     target = "French",
#'     candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija",
#'                    "Altai", "Denisova", "Chimp"),
#'     minimize = TRUE,
#'     nsources = 2,
#'     ncores = 2,
#'     fulloutput = FALSE
#' )
#'
#' # filter out models which can clearly be rejected
#' fits <- qpAdm_filter(models, p = 0.05)
#' }
#'
#' @export
qpAdm_filter <- function(x, p = 0.05) {
    check_type(x, "admixr_result")
    if (attr(x, "command") != "qpAdm_rotation") {
        stop("Filtering implemented only for results of the qpAdm rotation procedure",
             call. = FALSE)
    }

    if (length(x) == 3)
        proportions <- x$proportions
    else
        proportions <- x

    ## get positions of columns with estimated admixture proportions
    prop_columns <- stringr::str_which(names(proportions), "prop")

    ## find out rows/models for which all proportions are in [0, 1] and
    ## pvalue is larger than the required cutoff
    pvalue <- proportions$pvalue > p
    constr_0 <- proportions[, prop_columns] >= 0
    constr_1 <- proportions[, prop_columns] <= 1
    constr <- apply(pvalue & constr_0 & constr_1, 1, all)

    ## filter all three sub-tables to only those models that fit the criteria
    proportions <- dplyr::arrange(proportions[constr, ], -pvalue)

    if (length(x) == 3) {
        x$proportions <- proportions
        x$ranks <- x$ranks[x$ranks$model %in% proportions$model, ]
        x$subsets <- x$subsets[x$subsets$model %in% proportions$model, ]
        ## filter also only to relevant remaining log output information
        attr(x, "log_output") <- attr(x, "log_output")[unique(proportions$model)]
        return(x)
    } else {
        attr(proportions, "log_output") <- attr(x, "log_output")[unique(proportions$model)]
        return(proportions)
    }
}
