#' Sample embeddings at points
#'
#' `tessera_sample()` adds the embedding at each point to the points.
#' `tessera_probe()` also adds each point's status, which says why a point
#' has no embedding.
#'
#' Each point is read from the nearest stored pixel of its own UTM zone,
#' without interpolation. Near a zone seam, a point without an embedding is
#' also looked up in the neighbouring zone. Points are projected directly
#' from their own CRS to the zone's, and points already in the zone's CRS
#' are used as given.
#'
#' Points are read in batches of nearby points, of at most the option
#' `geotessera.batch_points` (default 10000), which bounds memory use. With
#' [mirai](https://mirai.r-lib.org) daemons set by [mirai::daemons()], batches
#' are read in parallel; each daemon opens the store itself and needs
#' geotessera installed. Most of the time goes to downloading, about 110 KB
#' for each point far from others, so a few daemons usually saturate the
#' connection. Give [tessera_open()] a `cache_dir` so that daemons share
#' downloads.
#'
#' @param x A store from [tessera_open()].
#' @param points The points to sample: a data frame with `lon` and `lat`
#'   (or `longitude` and `latitude`, or `x` and `y`) columns of WGS84
#'   coordinates, an `sf`, `sfc` or `SpatVector` point layer in any CRS, a
#'   two-column matrix of longitudes and latitudes, or a numeric vector
#'   `c(lon, lat)` for one point.
#' @param year The embedding year, one of `x$years`.
#' @param depth The number of leading bands of a published embedding
#'   prefix ("Matryoshka" depth), one of `names(x$depths)`, such as 16 on
#'   a v2 store. The default is the full embedding.
#' @param bands One-based band indices to return. The default is all bands
#'   at the selected depth.
#' @param cross_zone Whether to look up points near a seam in the
#'   neighbouring zone.
#' @param search_px The radius in pixels within which an unwritten pixel
#'   falls back to the nearest valid one. Zero disables the fallback. Open
#'   water is never searched past.
#' @return `points` with one numeric column per band added, named `band_1`,
#'   `band_2` and so on after the band indices, which are `NA` for points without an embedding. A
#'   data frame, tibble, `sf` or `SpatVector` keeps its class and columns;
#'   an `sfc` becomes an `sf` object; a matrix or vector becomes a data frame
#'   with `lon` and `lat` columns.
#'
#'   `tessera_probe()` also adds `status`, a factor with levels `"valid"`,
#'   `"water"` for open water, `"nodata"` for an unwritten pixel and
#'   `"outside"` beyond coverage.
#' @family readers
#' @examples
#' gt <- tessera_open(tessera_example())
#' pts <- data.frame(
#'   site = c("field", "lake", "outside"),
#'   lon = c(0.0085, 0.0069, 10),
#'   lat = c(52.0045, 52.0072, 52)
#' )
#' tessera_sample(gt, pts, 2024)
#'
#' probed <- tessera_probe(gt, pts, 2024)
#' table(probed$status)
#'
#' # A feature matrix for modelling.
#' emb <- as.matrix(tessera_sample(gt, pts, 2024)[paste0("band_", 1:4)])
#' @export
tessera_sample <- function(x, points, year, depth = NULL, bands = NULL,
                           cross_zone = TRUE, search_px = 1L) {
  res <- probe_points(x, points, year, depth, bands, cross_zone, search_px)
  add_columns(points, res$lonlat, as.data.frame(res$embeddings))
}

#' @rdname tessera_sample
#' @export
tessera_probe <- function(x, points, year, depth = NULL, bands = NULL,
                          cross_zone = TRUE, search_px = 1L) {
  res <- probe_points(x, points, year, depth, bands, cross_zone, search_px)
  cols <- as.data.frame(res$embeddings)
  cols$status <- factor(res$status, levels = c(VALID, WATER, NODATA, OUTSIDE))
  add_columns(points, res$lonlat, cols)
}

# Embeddings as a matrix with one row per point, and statuses.
probe_points <- function(x, points, year, depth, bands, cross_zone, search_px) {
  check_tessera(x)
  array <- embeddings_array(x, depth)
  bands <- check_bands(bands, array$bands)
  if (!isTRUE(cross_zone) && !isFALSE(cross_zone)) {
    abort("`cross_zone` must be TRUE or FALSE")
  }
  search_px <- check_count(search_px, "search_px", zero = TRUE)
  pts <- as_points(points)
  n <- nrow(pts$xy)
  values <- matrix(NA_real_, n, array$bands,
                   dimnames = list(NULL, paste0("band_", seq_len(array$bands))))
  status <- rep(OUTSIDE, n)

  batches <- point_batches(pts$lonlat, getOption("geotessera.batch_points", 10000L))
  parts <- lapply(batches, function(i) {
    list(xy = pts$xy[i, , drop = FALSE], lonlat = pts$lonlat[i, , drop = FALSE])
  })
  results <- if (length(parts) > 1L && use_mirai()) {
    m <- mirai::mirai_map(parts, probe_batch_worker, .args = list(
      url = x$url, cache_dir = x$cache_dir, crs = pts$crs, year = year,
      array = array, cross_zone = cross_zone, search_px = search_px
    ))
    if (interactive()) m[mirai::.progress, mirai::.stop] else m[mirai::.stop]
  } else {
    lapply(parts, function(p) {
      probe_batch(x, p, pts$crs, year, array, cross_zone, search_px)
    })
  }
  for (b in seq_along(batches)) {
    values[batches[[b]], ] <- results[[b]]$values
    status[batches[[b]]] <- results[[b]]$status
  }

  if (!is.null(bands)) values <- values[, bands, drop = FALSE]
  list(embeddings = values, status = status, lonlat = pts$lonlat)
}

# Row indices of points in batches of at most `size`, ordered by zone and
# then by half-degree cell, so that a batch's points share shards. Batches
# bound the memory of a read, which holds every chunk it needs at once.
point_batches <- function(lonlat, size) {
  n <- nrow(lonlat)
  size <- check_count(size, "geotessera.batch_points")
  if (n <= size) return(list(seq_len(n)))
  lon <- lonlat[, 1L]
  lat <- lonlat[, 2L]
  ord <- order(floor((lon + 180) / 6), floor(lat * 2), floor(lon * 2))
  unname(split(ord, ceiling(seq_len(n) / size)))
}

# Whether the user has set mirai daemons to sample batches in parallel.
use_mirai <- function() {
  requireNamespace("mirai", quietly = TRUE) && mirai::daemons_set()
}

# probe_batch() on a mirai daemon, which opens the store once and reuses it
# for later batches.
probe_batch_worker <- function(part, url, cache_dir, crs, year, array,
                               cross_zone, search_px) {
  key <- paste(url, cache_dir %||% "", sep = "\n")
  gt <- worker_stores[[key]]
  if (is.null(gt)) {
    gt <- tessera_open(url, cache_dir)
    assign(key, gt, envir = worker_stores)
  }
  probe_batch(gt, part, crs, year, array, cross_zone, search_px)
}

worker_stores <- new.env(parent = emptyenv())

# Embeddings and statuses of a batch of points, `part$xy` in `crs` with
# WGS84 coordinates `part$lonlat`.
probe_batch <- function(gt, part, crs, year, array, cross_zone, search_px) {
  lon <- part$lonlat[, 1L]
  lat <- part$lonlat[, 2L]
  n <- length(lon)
  values <- matrix(NA_real_, n, array$bands)
  status <- rep(OUTSIDE, n)
  seen_water <- seen_nodata <- rep(FALSE, n)

  finite <- is.finite(lon) & is.finite(lat)
  candidates <- matrix(NA_integer_, n, 3L)
  candidates[finite, 1L] <- zone_for_lon(lon[finite])
  if (cross_zone && any(finite)) {
    seams <- seam_neighbours(lon[finite])
    candidates[finite, 2L] <- seams[, "west"]
    candidates[finite, 3L] <- seams[, "east"]
  }

  for (pass in 1:3) {
    pending <- which(status != VALID & !is.na(candidates[, pass]))
    for (z in unique(candidates[pending, pass])) {
      sel <- pending[candidates[pending, pass] == z]
      zone <- zone_open(gt, z)
      if (is.null(zone)) next
      # Straight from the points' CRS, untouched when it is the zone's.
      en <- project_xy(part$xy[sel, 1L], part$xy[sel, 2L], crs, zone$crs)
      res <- probe_zone(gt, zone, en$x, en$y, year, search_px, array$name, array$bands)
      hit <- res$status == VALID
      values[sel[hit], ] <- res$values[hit, ]
      status[sel[hit]] <- VALID
      seen_water[sel] <- seen_water[sel] | res$status == WATER
      seen_nodata[sel] <- seen_nodata[sel] | res$status == NODATA
    }
  }

  missing <- status != VALID
  status[missing & seen_nodata] <- NODATA
  status[missing & seen_water] <- WATER
  list(values = values, status = status)
}

# `points` with the columns of `cols` added, keeping its class where it has
# one.
add_columns <- function(points, lonlat, cols) {
  if (inherits(points, "sfc")) {
    return(sf::st_sf(cols, geometry = points))
  }
  if (inherits(points, "SpatVector")) {
    clash <- intersect(names(cols), names(points))
    if (length(clash)) abort_clash(clash)
    return(cbind(points, cols))
  }
  if (!is.data.frame(points)) {
    return(cbind(data.frame(lon = lonlat[, 1L], lat = lonlat[, 2L]), cols))
  }
  clash <- intersect(names(cols), names(points))
  if (length(clash)) abort_clash(clash)
  points[names(cols)] <- cols
  points
}

abort_clash <- function(clash) {
  abort(sprintf("`points` already has columns named %s",
                paste0("`", clash, "`", collapse = ", ")))
}
