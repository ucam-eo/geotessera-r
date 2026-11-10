# Conditions and small helpers.

# Errors carry the class `geotessera_error`, so callers can catch them with
# tryCatch(..., geotessera_error = ).
abort <- function(..., class = NULL) {
  stop(errorCondition(paste0(...), class = c(class, "geotessera_error"), call = NULL))
}

warn <- function(..., class = NULL) {
  warning(warningCondition(paste0(...), class = c(class, "geotessera_warning"),
                           call = NULL))
}

check_installed <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    abort(sprintf("Package '%s' is required for %s", pkg, what))
  }
}

check_tessera <- function(x) {
  if (!inherits(x, "tessera")) {
    abort("`x` must be a store from tessera_open()")
  }
}
