#' Export a region as GeoTIFF files
#'
#' Writes one file per intersecting UTM zone, named
#' `tessera_YEAR_utmNN.tif`, in row strips. Files contain float32
#' embeddings on the native grid with NaN nodata, layer names `band_1`,
#' `band_2` and so on as returned by [tessera_region()], and `TESSERA_*`
#' metadata tags recording the year, source, model and bands. Each
#' completed file replaces its destination.
#'
#' Embeddings of different datasets must not be mixed. For a store listed
#' by [tessera_datasets()], each file is also tagged with the dataset
#' version and variant, and a `tessera_metadata.json` file in `output_dir`
#' records the dataset, as the Python `geotessera` package does. Exporting
#' into a directory that records another dataset is an error.
#'
#' @inheritParams tessera_region
#' @param bbox The bounds to export, in any form accepted by
#'   [tessera_region()]. Split bounds crossing the antimeridian into two
#'   requests.
#' @param output_dir The directory to write the files to.
#' @param strip_rows The maximum number of rows read at a time. Strips of
#'   wide regions are shorter, so that each holds about 64 MiB.
#' @param compress The GeoTIFF compression: one of `"lzw"`, `"deflate"`,
#'   `"zstd"` or `"none"`.
#' @param progress `TRUE` to show a text progress bar, or a function called
#'   as `progress(current, total, status)` after each strip is written.
#' @param dry_run Whether to return size estimates without reading
#'   embeddings or creating files.
#' @return The paths of the written files, invisibly. With
#'   `dry_run = TRUE`, a data frame with columns `zone`, `width`, `height` and
#'   `uncompressed_bytes`, which estimate the output before compression.
#' @family readers
#' @examples
#' gt <- tessera_open(tessera_example())
#' bbox <- c(-0.0067, 52.0018, 0.0087, 52.0036)
#' tessera_export_geotiffs(gt, bbox, 2024, tempdir(), dry_run = TRUE)
#'
#' files <- tessera_export_geotiffs(gt, bbox, 2024, tempdir(), bands = 1:3)
#' terra::rast(files[1])
#' unlink(files)
#' @export
tessera_export_geotiffs <- function(x, bbox, year, output_dir, depth = NULL,
                                    bands = NULL, strip_rows = 128L,
                                    compress = c("lzw", "deflate", "zstd", "none"),
                                    progress = NULL, dry_run = FALSE) {
  check_tessera(x)
  year <- check_year(x, year)
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir)) {
    abort("`output_dir` must be a single path")
  }
  compress <- match.arg(compress)
  strip_rows <- check_count(strip_rows, "strip_rows")
  array <- embeddings_array(x, depth)
  bands <- check_bands(bands, array$bands) %||% seq_len(array$bands)
  windows <- region_windows(x, as_bounds(bbox), year)
  if (!length(windows)) {
    abort("Region does not intersect any available zone grid")
  }
  estimates <- data.frame(
    zone = vapply(windows, function(w) w$zone$zone, integer(1)),
    width = vapply(windows, function(w) w$x1 - w$x0, numeric(1)),
    height = vapply(windows, function(w) w$y1 - w$y0, numeric(1))
  )
  estimates$uncompressed_bytes <- estimates$width * estimates$height * length(bands) * 4
  if (dry_run) return(estimates)

  dataset <- x$dataset
  if (!is.null(dataset)) check_dataset_dir(output_dir, dataset)
  for (w in windows) warn_incomplete(w$zone, year, w$y0, w$y1)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows_for <- function(w) strip_height(w$x1 - w$x0, length(bands), strip_rows)
  total <- sum(vapply(windows, function(w) ceiling((w$y1 - w$y0) / rows_for(w)),
                      numeric(1)))
  if (isTRUE(progress)) {
    bar <- utils::txtProgressBar(max = total, style = 3)
    on.exit(close(bar), add = TRUE)
    progress <- function(current, total, status) utils::setTxtProgressBar(bar, current)
  }
  completed <- 0
  files <- character(0)
  for (w in windows) {
    zone <- w$zone
    rows <- rows_for(w)
    path <- file.path(output_dir, sprintf("tessera_%d_utm%02d.tif", as.integer(year),
                                          zone$zone))
    tmp <- tempfile(".tessera-", tmpdir = output_dir, fileext = ".tif")
    r <- terra::rast(nrows = w$y1 - w$y0, ncols = w$x1 - w$x0, nlyrs = length(bands),
                     extent = window_extent(zone, w$x0, w$x1, w$y0, w$y1),
                     crs = zone$crs)
    names(r) <- paste0("band_", bands)
    terra::metags(r) <- c(
      if (!is.null(dataset)) dataset_tags(dataset),
      TESSERA_YEAR = as.character(year),
      TESSERA_SOURCE = x$url,
      TESSERA_MODEL = x$model_version,
      TESSERA_BUILD_VERSION = x$build_version,
      TESSERA_BANDS = paste(bands, collapse = ",")
    )
    ok <- FALSE
    tryCatch({
      terra::writeStart(r, tmp, overwrite = TRUE, datatype = "FLT4S", NAflag = NaN,
                        gdal = c(paste0("COMPRESS=", toupper(compress)), "TILED=YES",
                                 "BLOCKXSIZE=256", "BLOCKYSIZE=256", "BIGTIFF=IF_SAFER"))
      for (top in seq(w$y0, w$y1 - 1, by = rows)) {
        end <- min(top + rows, w$y1)
        q <- read_window(x, zone, array$name, w$ti, w$x0, w$x1, top, end, bands)
        values <- dequantize_array(q$embeddings, q$scales)
        # An [x, y, band] array flattens to row-major cells, one column per band.
        values <- matrix(values, ncol = length(bands))
        terra::writeValues(r, values, top - w$y0 + 1, end - top)
        completed <- completed + 1
        if (is.function(progress)) {
          progress(completed, total, sprintf("Writing UTM zone %d", zone$zone))
        }
      }
      terra::writeStop(r)
      if (file.exists(path)) unlink(path)
      if (!file.rename(tmp, path)) abort("Cannot replace ", path)
      ok <- TRUE
    }, finally = if (!ok) unlink(tmp))
    files <- c(files, path)
  }
  if (!is.null(dataset)) {
    write_dataset_metadata(output_dir, dataset,
                           list(format = "tiff", year = as.integer(year), source = x$url))
  }
  invisible(files)
}

# Native windows of every zone intersecting the bounds.
region_windows <- function(gt, bounds, year) {
  bbox <- bounds$lonlat
  west <- bbox[1L]; south <- bbox[2L]; east <- bbox[3L]; north <- bbox[4L]
  # An east edge on a zone boundary does not reach into the next zone.
  last <- as.integer(min(60, ceiling((east + 180) / 6)))
  out <- list()
  for (z in seq(zone_for_lon(west), last)) {
    zone <- zone_open(gt, z)
    if (is.null(zone)) next
    ti <- zone_time_index(zone, year)
    zone_bbox <- c(max(west, z * 6 - 186), south, min(east, z * 6 - 180), north)
    env <- utm_envelope(zone_bbox, zone$crs)
    if (!same_crs(bounds$crs, wgs84)) {
      # Projected bounds are not widened by the detour through WGS84.
      exact <- bounds_envelope(bounds, zone$crs)
      env <- c(max(env[1L], exact[1L]), max(env[2L], exact[2L]),
               min(env[3L], exact[3L]), min(env[4L], exact[4L]))
      if (env[1L] >= env[3L] || env[2L] >= env[4L]) next
    }
    win <- zone_window(zone, env)
    if (is.null(win)) next
    out[[length(out) + 1L]] <- c(win, list(zone = zone, ti = ti))
  }
  out
}
