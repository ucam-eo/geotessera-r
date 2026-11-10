# Per-store cache subdirectory. Published stores key by dataset path; other
# locations key by a slug and a digest of the location.
store_cache_key <- function(location) {
  location <- sub("/+$", "", location)
  prefix <- paste0(tessera_mirror_url, "/zarr/")
  if (startsWith(location, prefix)) {
    dataset <- substring(location, nchar(prefix) + 1L)
    if (nzchar(dataset) && !grepl("/", dataset, fixed = TRUE)) {
      return(gsub("[^A-Za-z0-9._-]+", "_", dataset))
    }
  }
  digest <- substr(digest::digest(location, "sha256", serialize = FALSE), 1, 12)
  slug <- gsub("[^A-Za-z0-9._-]+", "_", sub("^.*://", "", location))
  slug <- gsub("^_+|_+$", "", slug)
  paste0(substring(slug, max(1L, nchar(slug) - 79L)), "-", digest)
}

# A read-only HTTP store over Zarr v3 consolidated metadata that retries
# failed requests and fetches byte ranges concurrently. The zarr store
# returns the root metadata for every node of a v3 hierarchy. The root
# metadata persists in `cache_file`, and other whole objects, such as the
# small unsharded coordinate arrays, in the chunk cache.
tessera_http_store <- R6::R6Class(
  "tessera_http_store",
  inherit = zarr::zarr_httpstore,
  cloneable = FALSE,
  private = list(
    .cache_file = NULL,
    .cache = NULL,

    request = function(key, byte_range = NULL) {
      root_json <- is.null(byte_range) && identical(key, "zarr.json")
      cached <- private$.cache_file
      if (root_json && !is.null(cached) && file.exists(cached)) {
        return(readBin(cached, "raw", file.size(cached)))
      }
      object <- is.null(byte_range) && !root_json && !is.null(private$.cache)
      if (object) {
        body <- cache_get(private$.cache, paste0("object:", key))
        if (!is.null(body)) return(body)
      }
      body <- self$get_ranges(key, list(byte_range))[[1L]]
      if (object && !is.null(body)) cache_put(private$.cache, paste0("object:", key), body)
      if (root_json && !is.null(cached) && !is.null(body)) {
        dir.create(dirname(cached), recursive = TRUE, showWarnings = FALSE)
        tmp <- tempfile(tmpdir = dirname(cached))
        writeBin(body, tmp)
        file.rename(tmp, cached)
      }
      body
    }
  ),
  public = list(
    initialize = function(url, cache_file = NULL, cache = NULL) {
      private$.cache_file <- cache_file
      private$.cache <- cache
      super$initialize(url)
    },

    get_metadata = function(prefix) {
      path <- gsub("^/+|/+$", "", prefix)
      if (!nzchar(path)) private$.metadata
      else private$.metadata$consolidated_metadata$metadata[[path]]
    },

    get_ranges = function(keys, ranges) {
      http_get_ranges(paste(private$.base_url, keys, sep = "/"), ranges)
    }
  )
)

#' Open a Tessera Zarr store
#'
#' Opens a local or remote Zarr v3 store of Tessera embeddings and reads its
#' consolidated metadata. The returned handle is passed to every other
#' function and caches shard indices and compressed chunks in memory, so
#' reuse one handle for many queries.
#'
#' @param store A store URL from [tessera_store_url()] or the path of a
#'   local store. The default is the public v1.1 dClimate store.
#' @param cache_dir A directory in which to keep a remote store's metadata
#'   and downloaded chunks between sessions, or `NULL` to keep them only in
#'   memory. Each store URL uses its own subdirectory, whose chunks are
#'   bounded by the option `geotessera.disk_cache_bytes`. Cached data is
#'   never revalidated: delete a store's subdirectory to see changes to the
#'   store. The default is the option `geotessera.cache_dir`. A good choice
#'   is `tools::R_user_dir("geotessera", "cache")`.
#' @return A `tessera` object. It is an environment whose fields include:
#'   \describe{
#'     \item{`url`}{The store location.}
#'     \item{`years`}{The embedding years available, as integers.}
#'     \item{`n_bands`}{The number of bands in the full embedding.}
#'     \item{`depths`}{The published embedding depths, as the names of a
#'       character vector.}
#'     \item{`zone_ids`}{The UTM zones in the store.}
#'     \item{`model_version`, `build_version`}{Provenance of the
#'       embeddings.}
#'     \item{`dataset`}{The published dataset, as a list with `version`
#'       and `variant`, or `NULL` for a store not listed by
#'       [tessera_datasets()].}
#'   }
#' @family stores
#' @examples
#' gt <- tessera_open(tessera_example())
#' gt
#' gt$years
#'
#' \donttest{
#' # The public v1.1 store; requires an internet connection.
#' gt <- tessera_open()
#' gt$dataset
#' }
#' @export
tessera_open <- function(store = tessera_store_url(),
                         cache_dir = getOption("geotessera.cache_dir")) {
  if (!is.character(store) || length(store) != 1L || is.na(store)) {
    abort("`store` must be a URL or a path")
  }
  url <- sub("/+$", "", store)
  remote <- grepl("^https?://", url)
  if (!grepl("://", url, fixed = TRUE)) url <- normalizePath(url, mustWork = FALSE)
  if (!is.null(cache_dir)) cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
  store_dir <- if (remote && !is.null(cache_dir)) {
    file.path(cache_dir, store_cache_key(url))
  }
  cache <- new_cache(getOption("geotessera.cache_bytes", 256 * 1024^2),
                     if (!is.null(store_dir)) file.path(store_dir, "chunks"))
  zstore <- tryCatch(
    if (remote) {
      cache_file <- if (!is.null(store_dir)) file.path(store_dir, "zarr.json")
      tessera_http_store$new(url, cache_file = cache_file, cache = cache)
    } else if (grepl("://", url, fixed = TRUE)) {
      abort("Only http(s) URLs and local paths are supported")
    } else {
      zarr::zarr_localstore$new(url, read_only = TRUE)
    },
    error = function(e) {
      abort(sprintf(
        "Failed to open Zarr store at '%s': %s\nSee tessera_datasets() for the published stores.",
        url, conditionMessage(e)
      ), class = intersect(class(e), "geotessera_http_error"))
    }
  )

  root_meta <- zstore$get_metadata("/")
  if (is.null(root_meta) || !identical(root_meta$node_type, "group")) {
    abort(sprintf("No Zarr group at '%s'", url))
  }
  consolidated <- root_meta$consolidated_metadata$metadata
  if (!length(consolidated)) {
    abort(sprintf(
      "Zarr store at '%s' has no consolidated metadata; consolidate it with zarr.consolidate_metadata() in Python",
      url
    ))
  }
  root_meta$consolidated_metadata <- NULL
  attrs <- root_meta$attributes

  gt <- new.env(parent = emptyenv())
  gt$url <- url
  gt$cache_dir <- cache_dir
  gt$dataset <- dataset_for_location(url)
  gt$store <- zstore
  gt$consolidated <- consolidated
  gt$root <- zarr::zarr$new(zstore)$root
  gt$attributes <- attrs
  gt$model_version <- attrs[["geoemb:model"]] %||% ""
  gt$build_version <- attrs[["geoemb:build_version"]] %||% ""
  gt$n_bands <- as.integer(attrs[["geoemb:dimensions"]] %||% 128L)
  depths <- attrs[["geoemb:depths"]]
  gt$depths <- if (length(depths)) {
    stats::setNames(
      vapply(depths, function(d) as.character(d$array), character(1)),
      vapply(depths, function(d) as.character(d$dimensions), character(1))
    )
  } else {
    stats::setNames("embeddings", gt$n_bands)
  }
  gt$zones <- new.env(parent = emptyenv())
  gt$cache <- cache
  class(gt) <- "tessera"

  gt$zone_ids <- store_zone_ids(gt)
  if (!length(gt$zone_ids)) {
    abort(sprintf("No UTM zone groups (utmNN) in the Zarr store at '%s'", url))
  }
  gt$years <- integer(0)
  for (z in gt$zone_ids) {
    zone <- tryCatch(zone_open(gt, z), error = function(e) NULL)
    if (!is.null(zone)) {
      gt$years <- zone$years
      break
    }
  }
  gt
}

store_zone_ids <- function(gt) {
  ids <- grep("^utm[0-9]{2}$", names(gt$consolidated), value = TRUE)
  sort(as.integer(substring(ids, 4L)))
}

node_metadata <- function(gt, path) {
  gt$consolidated[[path]]
}

#' @export
format.tessera <- function(x, ...) {
  c(
    paste0("<tessera> ", x$url),
    if (!is.null(x$dataset)) {
      paste0("Dataset: v", x$dataset$version, " ", x$dataset$variant)
    },
    if (nzchar(x$model_version)) paste0("Model  : ", x$model_version),
    paste0("Years  : ", paste(x$years, collapse = ", ")),
    paste0("Zones  : ", length(x$zone_ids)),
    paste0("Depths : ", paste(names(x$depths), collapse = ", "))
  )
}

#' @export
print.tessera <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

# Array name and band count for a depth; NULL selects the full embedding.
embeddings_array <- function(gt, depth) {
  if (is.null(depth)) return(list(name = "embeddings", bands = gt$n_bands))
  key <- as.character(depth)
  if (length(depth) != 1L || !key %in% names(gt$depths)) {
    abort(sprintf("depth %s not in this store; available: %s",
                 paste(depth, collapse = ","),
                 paste(sort(as.integer(names(gt$depths))), collapse = ", ")))
  }
  list(name = gt$depths[[key]], bands = as.integer(depth))
}

# `year` as an integer, which must be one of the store's years.
check_year <- function(gt, year) {
  if (!is.numeric(year) || length(year) != 1L || !is.finite(year) ||
      year != round(year)) {
    abort("`year` must be a single whole number")
  }
  if (!year %in% gt$years) {
    abort(sprintf("Year %d is unavailable; available: %s", as.integer(year),
                  paste(gt$years, collapse = ", ")))
  }
  as.integer(year)
}

#' Path of the example store
#'
#' The package ships a small synthetic store for examples and tests. It
#' covers 256 by 256 pixels of 10 m either side of the seam between UTM zones
#' 30 and 31 at latitude 52 degrees north, for 2023 and 2024, with four
#' bands and a two-band depth. Its values are synthetic, not embeddings.
#'
#' @return The path of the example store, for [tessera_open()].
#' @family stores
#' @examples
#' gt <- tessera_open(tessera_example())
#' tessera_sample(gt, c(0.0085, 52.0045), 2024)
#' @export
tessera_example <- function() {
  system.file("extdata", "example.zarr", package = "geotessera", mustWork = TRUE)
}
