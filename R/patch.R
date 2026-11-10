#' Read a patch centred on a point
#'
#' Reads a square of `size_px` pixels whose pixel `[size_px %/% 2 + 1,
#' size_px %/% 2 + 1]` contains the point.
#'
#' The patch is laid on the native grid of one UTM zone, so its pixels are
#' copied from that zone without resampling: the zone whose CRS the point
#' is given in, or else the zone holding the point. When the patch crosses
#' into a neighbouring zone, the pixels whose longitude lies there are taken
#' from that zone by nearest-neighbour resampling, and each zone fills the
#' other's gaps.
#'
#' Setting `dst_crs` lays the patch on a grid in that CRS centred on the
#' point instead. If `dst_crs` is the CRS of one of the zones the patch
#' covers, the grid is that zone's native grid and only the other zones are
#' resampled; otherwise every pixel is resampled.
#'
#' @inheritParams tessera_region
#' @param point The patch centre, as one point in any form accepted by
#'   `points` in [tessera_sample()]: a numeric vector `c(lon, lat)`, a
#'   one-row data frame or matrix, or a one-point `sf`, `sfc` or
#'   `SpatVector` layer in any CRS.
#' @param size_px The patch width and height in pixels.
#' @param dst_crs A projected CRS with metre units for the output grid, or
#'   `NULL` for the native grid of the point's zone.
#' @param method The resampling method of [terra::project()] for pixels
#'   taken from a zone other than the one whose grid the patch uses. Only
#'   `"near"` leaves embeddings unblended.
#' @return A `SpatRaster` of `size_px` rows and columns, with `NA` where the
#'   store holds no embedding.
#' @family readers
#' @examples
#' gt <- tessera_open(tessera_example())
#' p <- tessera_patch(gt, c(0.0175, 52.0011), 2024, size_px = 16)
#' p
#'
#' # A patch across the seam between zones 30 and 31, on zone 31's grid.
#' seam <- tessera_patch(gt, c(0, 52), 2024, size_px = 24)
#' terra::crs(seam, describe = TRUE)$code
#'
#' # The same patch on zone 30's grid, from a point in its CRS.
#' pt <- terra::project(terra::vect(cbind(0, 52), crs = "EPSG:4326"), "EPSG:32630")
#' terra::crs(tessera_patch(gt, pt, 2024, size_px = 24), describe = TRUE)$code
#' @export
tessera_patch <- function(x, point, year, size_px, depth = NULL, bands = NULL,
                          dst_crs = NULL, method = "near") {
  check_tessera(x)
  year <- check_year(x, year)
  pt <- as_points(point)
  if (nrow(pt$xy) != 1L) abort("`point` must be a single point")
  lon <- pt$lonlat[1L, 1L]
  lat <- pt$lonlat[1L, 2L]
  if (!is.finite(lon) || !is.finite(lat) || abs(lon) > 180 || abs(lat) > 90) {
    abort("`point` must lie within WGS84 bounds")
  }
  size_px <- check_count(size_px, "size_px")
  if (!is.null(dst_crs)) {
    dst_crs <- terra::crs(dst_crs)
    if (!nzchar(dst_crs) || terra::is.lonlat(dst_crs) ||
        !isTRUE(all.equal(terra::linearUnits(terra::rast(crs = dst_crs)), 1))) {
      abort("`dst_crs` must be a projected CRS with metre units")
    }
  }
  array <- embeddings_array(x, depth)
  bands <- check_bands(bands, array$bands)
  band_ids <- bands %||% seq_len(array$bands)

  centre <- zone_open(x, zone_for_lon(lon))
  if (is.null(centre)) {
    abort(sprintf("Zone %d is not in this store", zone_for_lon(lon)))
  }
  px <- centre$px
  half <- size_px * px / 2
  c0 <- project_xy(lon, lat, wgs84, centre$crs)
  corners <- project_xy(c0$x + c(-half, -half, half, half),
                        c0$y + c(-half, half, -half, half), centre$crs, wgs84)
  zones <- zones_spanned(corners$x, lon)

  # The zone whose native grid the patch uses, if any: the zone named by
  # `dst_crs`, else the zone whose CRS the point is in, else the point's.
  grid <- if (is.null(dst_crs)) centre
  want <- dst_crs %||% pt$crs
  for (z in zones) {
    zone <- zone_open(x, z)
    if (!is.null(zone) && same_crs(want, zone$crs)) {
      grid <- zone
      break
    }
  }
  if (!is.null(grid)) {
    # Straight from the point's CRS, untouched when it is the grid's.
    ce <- project_xy(pt$xy[1L, 1L], pt$xy[1L, 2L], pt$crs, grid$crs)
    win <- patch_window(grid, ce$x, ce$y, size_px)
    if (length(zones) == 1L) {
      return(patch_native(x, grid, win$x0, win$y0, year, size_px, array, bands, band_ids))
    }
    target <- grid$crs
    extent <- window_extent(grid, win$x0, win$x0 + size_px, win$y0, win$y0 + size_px)
  } else {
    target <- dst_crs
    ce <- project_xy(pt$xy[1L, 1L], pt$xy[1L, 2L], pt$crs, target)
    ox <- ce$x - (size_px %/% 2 + 0.5) * px
    oy <- ce$y + (size_px %/% 2 + 0.5) * px
    extent <- terra::ext(ox, ox + size_px * px, oy - size_px * px, oy)
  }
  patch_merged(x, zones, target, extent, year, size_px, array, bands, band_ids, method)
}

# The 0-based pixel offsets of a patch whose centre pixel holds (ce, cn).
patch_window <- function(zone, ce, cn, size_px) {
  ix <- ceiling((ce - zone_x(zone, 0)) / zone$px - 0.5)
  iy <- ceiling((zone_y(zone, 0) - cn) / zone$px - 0.5)
  list(x0 = ix - size_px %/% 2, y0 = iy - size_px %/% 2)
}

# The patch at pixel offsets (x0, y0) of a zone grid, NA beyond it.
patch_native <- function(gt, zone, x0, y0, year, size_px, array, bands, band_ids) {
  ti <- zone_time_index(zone, year)
  cx0 <- max(0, x0)
  cx1 <- min(zone$width, x0 + size_px)
  cy0 <- max(0, y0)
  cy1 <- min(zone$height, y0 + size_px)
  out <- array(NA_real_, c(size_px, size_px, length(band_ids)))
  if (cx1 > cx0 && cy1 > cy0) {
    q <- read_window(gt, zone, array$name, ti, cx0, cx1, cy0, cy1, bands)
    out[(cx0 - x0 + 1):(cx1 - x0), (cy0 - y0 + 1):(cy1 - y0), ] <-
      dequantize_array(q$embeddings, q$scales)
  }
  as_spatraster(out, window_extent(zone, x0, x0 + size_px, y0, y0 + size_px),
                zone$crs, band_ids)
}

# Pixel offsets of an extent on a zone's grid in the same CRS, or NULL when
# the extent is not aligned with the grid.
grid_offsets <- function(zone, extent, target) {
  if (!same_crs(zone$crs, target)) return(NULL)
  dx <- (extent$xmin - zone$ox) / zone$px
  dy <- (zone$oy - extent$ymax) / zone$px
  if (abs(dx - round(dx)) > 1e-6 || abs(dy - round(dy)) > 1e-6) return(NULL)
  list(x0 = round(dx), y0 = round(dy))
}

patch_merged <- function(gt, zones, target, extent, year, size_px, array, bands,
                         band_ids, method) {
  template <- terra::rast(nrows = size_px, ncols = size_px, nlyrs = length(band_ids),
                          extent = extent, crs = target)
  ox <- extent$xmin
  oy <- extent$ymax
  px <- (extent$xmax - extent$xmin) / size_px
  side <- size_px * px

  # Densified outline: straight edges curve in a zone's CRS.
  steps <- seq(0, side, length.out = 33L)
  outline_x <- c(ox + steps, ox + steps, rep(ox, 33L), rep(ox + side, 33L))
  outline_y <- c(rep(oy, 33L), rep(oy - side, 33L), oy - steps, oy - steps)

  centres <- (seq_len(size_px) - 0.5) * px
  gx <- rep(ox + centres, each = size_px)
  gy <- rep(oy - centres, times = size_px)
  owner <- matrix(zone_for_lon(project_xy(gx, gy, target, wgs84)$x), size_px, size_px)

  nb <- length(band_ids)
  out <- array(NA_real_, c(size_px, size_px, nb))
  owned <- matrix(FALSE, size_px, size_px)
  spare <- array(NA_real_, c(size_px, size_px, nb))
  spared <- matrix(FALSE, size_px, size_px)
  for (z in zones) {
    zone <- zone_open(gt, z)
    if (is.null(zone)) next
    native <- grid_offsets(zone, extent, target)
    if (!is.null(native)) {
      # The patch lies on this zone's grid: copy its pixels unchanged.
      moved <- terra::as.array(patch_native(gt, zone, native$x0, native$y0, year,
                                            size_px, array, bands, band_ids))
    } else {
      ti <- zone_time_index(zone, year)
      outline <- project_xy(outline_x, outline_y, target, zone$crs)
      pad <- 2 * zone$px
      win <- zone_window(zone, c(min(outline$x) - pad, min(outline$y) - pad,
                                 max(outline$x) + pad, max(outline$y) + pad))
      if (is.null(win)) next
      q <- read_window(gt, zone, array$name, ti, win$x0, win$x1, win$y0, win$y1, bands)
      src <- as_spatraster(dequantize_array(q$embeddings, q$scales),
                           window_extent(zone, win$x0, win$x1, win$y0, win$y1),
                           zone$crs, band_ids)
      moved <- terra::as.array(terra::project(src, template, method = method))
    }
    valid <- rowSums(!is.na(moved), dims = 2L) > 0
    mine <- valid & owner == z
    out[rep(mine, nb)] <- moved[rep(mine, nb)]
    owned <- owned | mine
    extra <- valid & owner != z & !spared
    spare[rep(extra, nb)] <- moved[rep(extra, nb)]
    spared <- spared | extra
  }
  fill <- !owned & spared
  out[rep(fill, nb)] <- spare[rep(fill, nb)]

  r <- terra::rast(out, crs = target, extent = extent)
  names(r) <- paste0("band_", band_ids)
  r
}
