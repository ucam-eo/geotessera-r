# Published datasets, the Zarr stores serving them, and the provenance that
# keeps embeddings of different datasets apart.

tessera_mirror_url <- "https://data.source.coop/tessera/tessera"


# One inference run of a Tessera model per row. Embeddings of different
# rows do not interoperate, even within one version. The first row of a
# version is its default variant. `path` is the store under the mirror's
# zarr/ tree.
datasets <- data.frame(
  version = c("1.0", "1.1", "1.1", "2.0", "2.0"),
  variant = c("vultr", "dclimate", "cambridge", "2B-L~beta1", "2B-L~beta2"),
  description = c("Legacy model; frozen", "Complete global run", "Cambridge test run",
                  "v2 beta", "v2 beta"),
  path = c("v1", "v1.1-dclimate", "v1.1", "v2-2B-L~beta1", "v2-2B-L~beta2")
)

#' Published Tessera datasets
#'
#' Lists the datasets published on the public Tessera mirror. Each is one
#' inference run of a model version, a variant. Embeddings of different
#' datasets are unrelated, even within one version, so use one dataset per
#' analysis.
#'
#' @return A data frame with one row per dataset and columns `version`,
#'   `variant`, `description`, `default` (whether the variant is its
#'   version's default) and `url`, the store for [tessera_open()].
#' @family stores
#' @examples
#' tessera_datasets()
#' @export
tessera_datasets <- function() {
  data.frame(
    datasets[c("version", "variant", "description")],
    default = !duplicated(datasets$version),
    url = paste0(tessera_mirror_url, "/zarr/", datasets$path)
  )
}

#' URL of a published Tessera Zarr store
#'
#' Builds the URL of a dataset's store on the public Tessera mirror, for
#' [tessera_open()]. See [tessera_datasets()] for the published versions
#' and variants.
#'
#' `"1"`, `"1.0"` and `"v1.0"` all name v1. A string that is not a version,
#' such as `"v1.1-dclimate"`, is taken as a store path under the mirror.
#'
#' @param version A dataset version such as `"v1"`, `"v1.1"` or `"v2"`, or
#'   a store path under the Tessera mirror. The default is v1.1.
#' @param variant A variant of `version`, such as `"cambridge"`. The
#'   default is the version's default variant, `"dclimate"` for v1.1.
#' @return A URL, as a string.
#' @family stores
#' @examples
#' tessera_store_url()
#' tessera_store_url("v1.1", "cambridge")
#' tessera_store_url("v2")
#' @export
tessera_store_url <- function(version = "v1.1", variant = NULL) {
  if (!is.character(version) || length(version) != 1L || is.na(version)) {
    abort("`version` must be a single string")
  }
  if (!is.null(variant) &&
      (!is.character(variant) || length(variant) != 1L || is.na(variant))) {
    abort("`variant` must be a single string")
  }
  norm <- dataset_version_norm(version)
  if (is.na(norm)) {
    if (!is.null(variant)) abort("`variant` needs a version such as \"v1.1\"")
    return(paste0(tessera_mirror_url, "/zarr/", version))
  }
  known <- datasets[datasets$version == norm, ]
  if (!nrow(known)) {
    path <- dataset_version_path(norm)
    if (!is.null(variant)) path <- paste0(path, "-", variant)
    return(paste0(tessera_mirror_url, "/zarr/", path))
  }
  variant <- variant %||% known$variant[1L]
  row <- known[known$variant == variant, ]
  if (!nrow(row)) {
    abort(sprintf("v%s has no variant '%s'; available: %s", norm, variant,
                  paste(known$variant, collapse = ", ")))
  }
  paste0(tessera_mirror_url, "/zarr/", row$path)
}

# "1.0" for "v1", "1", "1.0" or "v1.0"; NA for strings that are not versions.
dataset_version_norm <- function(version) {
  m <- regmatches(version, regexec("^v?([0-9]+)(?:\\.([0-9]+))?$", version))[[1L]]
  if (!length(m)) return(NA_character_)
  paste0(m[2L], ".", if (nzchar(m[3L])) m[3L] else "0")
}

# The mirror's name for a normalised version: "v1" for "1.0", "v1.1" for
# "1.1".
dataset_version_path <- function(norm) {
  ifelse(endsWith(norm, ".0"), paste0("v", sub("\\.0$", "", norm)), paste0("v", norm))
}

# The published dataset served at a store location, as a list with
# `version` and `variant`, or NULL for any other store.
dataset_for_location <- function(location) {
  location <- sub("/+$", "", location)
  i <- match(location, paste0(tessera_mirror_url, "/zarr/", datasets$path))
  if (is.na(i)) return(NULL)
  list(version = datasets$version[i], variant = datasets$variant[i])
}

# GeoTIFF tags identifying a dataset, as written by the Python package.
dataset_tags <- function(dataset) {
  c(TESSERA_DATASET_VERSION = dataset$version,
    TESSERA_DATASET_VERSION_PATH = dataset_version_path(dataset$version),
    TESSERA_DATASET_VARIANT = dataset$variant)
}

# A sidecar in each output directory records the dataset that populated
# it, in the format of the Python package, so that exports of different
# datasets are never mixed in one directory.
metadata_file <- "tessera_metadata.json"

# The dataset recorded in a directory's sidecar, or NULL.
recorded_dataset <- function(dir) {
  path <- file.path(dir, metadata_file)
  if (!file.exists(path)) return(NULL)
  meta <- tryCatch(jsonlite::read_json(path), error = function(e) NULL)
  version <- meta$dataset_version_path %||% meta$dataset_version
  variant <- meta$dataset_variant
  if (!is.character(version) || !is.character(variant)) return(NULL)
  norm <- dataset_version_norm(version)
  if (is.na(norm)) return(NULL)
  list(version = norm, variant = variant)
}

check_dataset_dir <- function(dir, dataset) {
  recorded <- recorded_dataset(dir)
  if (!is.null(recorded) && !identical(recorded, dataset)) {
    abort(sprintf(
      "'%s' holds v%s %s embeddings, but this store is v%s %s. Embeddings of different datasets cannot be interchanged; use a separate directory.",
      dir, recorded$version, recorded$variant, dataset$version, dataset$variant
    ), class = "geotessera_dataset_error")
  }
}

write_dataset_metadata <- function(dir, dataset, extra = list()) {
  meta <- c(list(
    dataset_version = dataset$version,
    dataset_version_path = dataset_version_path(dataset$version),
    dataset_variant = dataset$variant,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00", tz = "UTC")
  ), extra)
  meta <- meta[order(names(meta))]
  path <- file.path(dir, metadata_file)
  tmp <- tempfile(".tessera-", tmpdir = dir, fileext = ".json")
  jsonlite::write_json(meta, tmp, auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    abort("Cannot write ", path)
  }
  invisible(path)
}
