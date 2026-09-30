#' Read a region
#'
#' `tessera_region()` reads the embeddings within bounds on the native grid
#' of one UTM zone. `tessera_region_quantized()` reads the stored
#' quantized values and per-pixel scales, a smaller representation that
#' [tessera_dequantize()] converts.
#'
#' Pixels are never resampled: the bounds select the stored pixels whose
#' centres they enclose. One zone serves the whole region, the zone whose
#' CRS the bounds are given in or else the zone holding their centre. Use
#' [tessera_patch()] to merge across zones, or [tessera_export_geotiffs()]
#' to export every intersecting zone.
#'
#' Bounds in the zone's own CRS are used as given, so a `SpatRaster` on the
#' zone grid, such as the result of an earlier call, selects exactly its own
#' cells. Bounds in any other CRS, including longitude and latitude, select
#' the extent of their projected outline, which can be wider than the bounds
#' themselves where the zone grid is rotated from true north.
#'
#' @inheritParams tessera_sample
#' @param bbox The bounds to read: a numeric vector of WGS84
#'   `c(west, south, east, north)`, a named vector with `xmin`, `ymin`,
#'   `xmax` and `ymax`, a WGS84 `SpatExtent`, or an `sf`, `sfc`, `bbox`,
#'   `SpatVector` or `SpatRaster` object in any CRS, whose extent is used.
#' @param bands One-based band indices to read. The default is all bands at
#'   the selected depth.
#' @return `tessera_region()` returns a `SpatRaster` in the zone's UTM CRS
#'   with one layer per band and `NA` where the store holds no embedding.
#'
#'   `tessera_region_quantized()` returns a `SpatRaster` on the same grid
#'   with the integer band layers followed by a `scale` layer, which is
#'   `NaN` for open water and `Inf` for an unwritten pixel.
#' @family readers
#' @examples
#' gt <- tessera_open(tessera_example())
#' bbox <- c(0.004, 52.005, 0.013, 52.009)
#' r <- tessera_region(gt, bbox, 2024)
#' r
#' terra::plotRGB(terra::stretch(r[[1:3]]))
#'
#' q <- tessera_region_quantized(gt, bbox, 2024)
#' names(q)
#' all.equal(tessera_dequantize(q), r)
#' @export
tessera_region <- function(x, bbox, year, depth = NULL, bands = NULL) {
  win <- region_setup(x, bbox, year, depth, bands)
  q <- read_window(x, win$zone, win$array, win$ti, win$x0, win$x1, win$y0, win$y1,
                   win$bands)
  as_spatraster(dequantize_array(q$embeddings, q$scales), win$extent, win$zone$crs,
                win$band_ids)
}

#' @rdname tessera_region
#' @export
tessera_region_quantized <- function(x, bbox, year, depth = NULL, bands = NULL) {
  win <- region_setup(x, bbox, year, depth, bands)
  q <- read_window(x, win$zone, win$array, win$ti, win$x0, win$x1, win$y0, win$y1,
                   win$bands)
  nb <- length(win$band_ids)
  values <- array(c(q$embeddings, q$scales), c(dim(q$scales), nb + 1L))
  r <- as_spatraster(values, win$extent, win$zone$crs, c(win$band_ids, NA))
  names(r) <- c(paste0("band_", win$band_ids), "scale")
  r
}

#' Apply a function to a region in row strips
#'
#' Reads the pixels of [tessera_region()] top to bottom in strips, so the
#' region need not fit in memory. Each strip is a `SpatRaster` on the same
#' grid as the whole region.
#'
#' @inheritParams tessera_region
#' @param FUN A function called with each strip as a `SpatRaster`,
#'   followed by `...`.
#' @param ... Further arguments to `FUN`.
#' @param strip_rows The maximum number of rows per strip. Strips of wide
#'   regions are shorter, so that each holds about 64 MiB.
#' @return A list of the values returned by `FUN`, one per strip.
#' @family readers
#' @examples
#' gt <- tessera_open(tessera_example())
#' bbox <- c(0.004, 52.005, 0.013, 52.009)
#' means <- tessera_map_region(gt, bbox, 2024, terra::global, "mean",
#'                             na.rm = TRUE, strip_rows = 16)
#' length(means)
#' @export
tessera_map_region <- function(x, bbox, year, FUN, ..., depth = NULL, bands = NULL,
                               strip_rows = 128L) {
  FUN <- match.fun(FUN)
  strip_rows <- check_count(strip_rows, "strip_rows")
  win <- region_setup(x, bbox, year, depth, bands)
  rows <- strip_height(win$x1 - win$x0, length(win$band_ids), strip_rows)
  tops <- seq(win$y0, win$y1 - 1, by = rows)
  lapply(tops, function(top) {
    end <- min(top + rows, win$y1)
    q <- read_window(x, win$zone, win$array, win$ti, win$x0, win$x1, top, end,
                     win$bands)
    strip <- as_spatraster(dequantize_array(q$embeddings, q$scales),
                           window_extent(win$zone, win$x0, win$x1, top, end),
                           win$zone$crs, win$band_ids)
    FUN(strip, ...)
  })
}

region_setup <- function(gt, bbox, year, depth, bands) {
  check_tessera(gt)
  bounds <- as_bounds(bbox)
  array <- embeddings_array(gt, depth)
  bands <- check_bands(bands, array$bands)
  lonlat <- bounds$lonlat
  centre <- (lonlat[1L] + lonlat[3L]) / 2
  spanned <- zones_spanned(c(lonlat[1L], lonlat[3L]), centre)
  # Bounds in the CRS of a zone they touch are read from that zone's grid.
  z <- zone_for_lon(centre)
  for (candidate in spanned) {
    zone <- zone_open(gt, candidate)
    if (!is.null(zone) && same_crs(bounds$crs, zone$crs)) {
      z <- candidate
      break
    }
  }
  if (length(spanned) > 1L) {
    warn(sprintf(
      "bbox spans UTM zones %s; reading zone %d only. Use tessera_patch() or tessera_export_geotiffs() for more than one zone.",
      paste(spanned, collapse = ", "), z))
  }
  zone <- zone_open(gt, z)
  if (is.null(zone)) abort(sprintf("Zone %d is not in this store", z))
  ti <- zone_time_index(zone, year)
  win <- zone_window(zone, bounds_envelope(bounds, zone$crs))
  if (is.null(win)) abort("Region does not intersect the zone grid")
  warn_incomplete(zone, year, win$y0, win$y1)
  c(win, list(
    zone = zone, ti = ti, array = array$name, bands = bands,
    band_ids = bands %||% seq_len(array$bands),
    extent = window_extent(zone, win$x0, win$x1, win$y0, win$y1)
  ))
}

check_bands <- function(bands, count) {
  if (is.null(bands)) return(NULL)
  if (!is.numeric(bands) || !length(bands) || any(!is.finite(bands)) ||
      any(bands != round(bands)) || any(bands < 1 | bands > count)) {
    abort(sprintf("`bands` must be one-based indices no greater than %d", count))
  }
  as.integer(bands)
}

check_count <- function(x, name, zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x) ||
      x < if (zero) 0 else 1) {
    abort(sprintf("`%s` must be a %s integer", name,
                  if (zero) "non-negative" else "positive"))
  }
  as.integer(x)
}

# Rows per strip of a window `width` pixels wide, at most `strip_rows`,
# bounding each strip's working memory, integer and double copies of every
# band, to about 64 MiB.
strip_height <- function(width, n_bands, strip_rows) {
  as.integer(min(strip_rows, max(1, (64 * 1024^2) %/% (width * (n_bands * 20 + 16)))))
}

# A SpatRaster from an [x, y, band] array.
as_spatraster <- function(values, extent, crs, band_ids) {
  r <- terra::rast(aperm(values, c(2L, 1L, 3L)), crs = crs, extent = extent)
  names(r) <- paste0("band_", band_ids)
  r
}
