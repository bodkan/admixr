# Create either specified or a temporary directory.
get_dir <- function(dir_name=NULL) {
    if (!is.null(dir_name)) {
        dir.create(dir_name)
    } else {
        dir_name <- tempdir()
    }

    dir_name
}

# Generate paths to the population file, parameter file and log file
# based on a specified directory.
get_files <- function(dir_name, prefix) {
    directory <- get_dir(dir_name)
    list(
        pop_file=file.path(directory, paste0(prefix, ".pop")),
        par_file=file.path(directory, paste0(prefix, ".par")),
        log_file=file.path(directory, paste0(prefix, ".log"))
    )
}
