# UTM zone groups: grid geometry, year lookup and array readers.

VALID <- "valid"
WATER <- "water"
NODATA <- "nodata"
OUTSIDE <- "outside"

# The zone group `utmNN`, or NULL when the store has no such zone.
zone_open <- function(gt, zone) {
  key <- as.character(zone)
  if (exists(key, envir = gt$zones, inherits = FALSE)) {
    z <- get(key, envir = gt$zones)
    return(if (isFALSE(z)) NULL else z)
  }
  path <- sprintf("utm%02d", as.integer(zone))
  meta <- node_metadata(gt, path)
  if (is.null(meta) || !identical(meta$node_type, "group")) {
    assign(key, FALSE, envir = gt$zones)
    return(NULL)
  }
  attrs <- meta$attributes
  transform <- as.numeric(attrs[["spatial:transform"]])
  code <- attrs[["proj:code"]]
  if (length(transform) != 6L || is.null(code)) {
    abort(sprintf("Zone '%s' lacks spatial:transform or proj:code", path))
  }
  group <- zarr::zarr_group$new(path, meta, gt$root, no_create_check = TRUE)

  scales <- node_metadata(gt, paste0(path, "/scales"))
  time_meta <- node_metadata(gt, paste0(path, "/time"))
  time <- zarr::zarr_array$new("time", time_meta, group, no_create_check = TRUE)
  time$raw_read <- TRUE

  z <- new.env(parent = emptyenv())
  z$zone <- as.integer(zone)
  z$path <- path
  z$group <- group
  z$crs <- code
  z$px <- transform[1L]
  z$ox <- transform[3L]
  z$oy <- transform[6L]
  z$height <- as.numeric(scales$shape[2L])
  z$width <- as.numeric(scales$shape[3L])
  z$years <- as.integer(time$read())
  z$source_nodata <- identical(attrs[["geotessera:mask_source"]], "source_nodata")
  z$incomplete <- incomplete_years(attrs[["geotessera:source_groups"]], z$years)
  z$readers <- new.env(parent = emptyenv())
  assign(key, z, envir = gt$zones)
  z
}

# Years each source group (such as 31N) holds but does not list as
# complete, by group. Incomplete years read as nodata.
incomplete_years <- function(groups, years) {
  out <- list()
  for (name in names(groups)) {
    complete <- groups[[name]]$years_complete
    if (is.null(complete)) next
    missing <- setdiff(years, as.integer(unlist(complete)))
    if (length(missing)) out[[name]] <- missing
  }
  out
}

# Warn when rows [y0, y1) of a zone reach a source group without a
# complete `year`. Groups ending N hold rows above the equator and groups
# ending S rows below it.
warn_incomplete <- function(zone, year, y0, y1) {
  top <- zone$oy - y0 * zone$px
  bottom <- zone$oy - y1 * zone$px
  for (name in names(zone$incomplete)) {
    reaches <- if (endsWith(name, "S")) bottom < 0 else top > 0
    if (reaches && as.integer(year) %in% zone$incomplete[[name]]) {
      warn(sprintf("%s has no complete %d; its pixels read as nodata", name,
                   as.integer(year)))
    }
  }
}

zone_reader <- function(gt, zone, name) {
  reader <- zone$readers[[name]]
  if (is.null(reader)) {
    meta <- node_metadata(gt, paste0(zone$path, "/", name))
    if (is.null(meta)) {
      abort(sprintf("Array '%s/%s' not in this store", zone$path, name))
    }
    arr <- zarr::zarr_array$new(name, meta, zone$group, no_create_check = TRUE)
    reader <- sharded_reader(arr, gt$store, gt$cache)
    assign(name, reader, envir = zone$readers)
  }
  reader
}

zone_time_index <- function(zone, year) {
  i <- match(as.integer(year), zone$years)
  if (length(year) != 1L || is.na(i)) {
    abort(sprintf("Year %s is unavailable; available: %s",
                 paste(year, collapse = ","), paste(zone$years, collapse = ", ")))
  }
  i - 1L
}

# Pixel centre coordinates for 0-based column and row indices.
zone_x <- function(zone, i) zone$ox + (i + 0.5) * zone$px
zone_y <- function(zone, j) zone$oy - (j + 0.5) * zone$px

# Nearest 0-based pixel index along an axis; ties choose the lower index.
nearest_index <- function(v, origin, step, n) {
  k <- floor((v - origin) / step - 0.5)
  best <- rep(NA_real_, length(v))
  dist <- rep(Inf, length(v))
  for (o in -1:2) {
    c <- pmin(pmax(k + o, 0), n - 1)
    d <- abs(origin + (c + 0.5) * step - v)
    better <- !is.na(d) & (d < dist | (d == dist & c < best))
    best[better] <- c[better]
    dist[better] <- d[better]
  }
  best
}

# Number of leading indices satisfying a monotone predicate.
count_leading <- function(pred, guess, n) {
  k <- min(max(guess, 0), n)
  while (k > 0 && !pred(k - 1)) k <- k - 1
  while (k < n && pred(k)) k <- k + 1
  k
}

# The pixel window whose centres lie within projected bounds, or NULL.
zone_window <- function(zone, ebox) {
  west <- min(ebox[1L], ebox[3L])
  east <- max(ebox[1L], ebox[3L])
  south <- min(ebox[2L], ebox[4L])
  north <- max(ebox[2L], ebox[4L])
  px <- zone$px
  x0 <- count_leading(function(i) zone_x(zone, i) < west,
                      ceiling((west - zone$ox) / px - 0.5), zone$width)
  x1 <- count_leading(function(i) zone_x(zone, i) <= east,
                      floor((east - zone$ox) / px - 0.5) + 1, zone$width)
  y0 <- count_leading(function(j) zone_y(zone, j) > north,
                      ceiling((zone$oy - north) / px - 0.5), zone$height)
  y1 <- count_leading(function(j) zone_y(zone, j) >= south,
                      floor((zone$oy - south) / px - 0.5) + 1, zone$height)
  if (x1 <= x0 || y1 <= y0) return(NULL)
  list(x0 = x0, x1 = x1, y0 = y0, y1 = y1)
}

window_extent <- function(zone, x0, x1, y0, y1) {
  terra::ext(zone$ox + x0 * zone$px, zone$ox + x1 * zone$px,
             zone$oy - y1 * zone$px, zone$oy - y0 * zone$px)
}

# Quantized embeddings [x, y, band] and scales [x, y] for a window.
read_window <- function(gt, zone, array, ti, x0, x1, y0, y1, bands = NULL) {
  emb_reader <- zone_reader(gt, zone, array)
  b0 <- 0
  b1 <- emb_reader$shape[2L]
  if (!is.null(bands)) {
    b0 <- min(bands) - 1
    b1 <- max(bands)
  }
  emb <- read_block(emb_reader, c(ti, b0, y0, x0), c(ti + 1, b1, y1, x1))
  dim(emb) <- dim(emb)[1:3]
  if (!is.null(bands) && !identical(as.numeric(bands), seq(b0 + 1, b1))) {
    emb <- emb[, , bands - b0, drop = FALSE]
  }
  scales <- read_block(zone_reader(gt, zone, "scales"), c(ti, y0, x0), c(ti + 1, y1, x1))
  dim(scales) <- dim(scales)[1:2]
  list(embeddings = emb, scales = scales)
}

#' Dequantize Tessera embeddings
#'
#' Tessera stores embeddings as 8-bit integers with a scale per pixel. This
#' multiplies each band by its pixel's scale, as [tessera_region()] does.
#'
#' @param x A `SpatRaster` from [tessera_region_quantized()].
#' @return A `SpatRaster` of the band layers of `x` multiplied by its
#'   `scale` layer. Pixels with a non-finite scale, which marks open water or
#'   an unwritten pixel, are `NA`.
#' @family readers
#' @examples
#' gt <- tessera_open(tessera_example())
#' q <- tessera_region_quantized(gt, c(0.004, 52.005, 0.013, 52.009), 2024)
#' tessera_dequantize(q)
#' @export
tessera_dequantize <- function(x) {
  if (!inherits(x, "SpatRaster") || !"scale" %in% names(x)) {
    abort("`x` must be a SpatRaster from tessera_region_quantized()")
  }
  bands <- names(x) != "scale"
  scale <- x[["scale"]]
  out <- x[[which(bands)]] * terra::ifel(is.finite(scale), scale, NA)
  names(out) <- names(x)[bands]
  out
}

# Dequantized values of an embeddings array whose leading dimensions are
# those of `scales`. Non-finite scales give NA.
dequantize_array <- function(embeddings, scales) {
  s <- as.vector(scales)
  s[!is.finite(s)] <- NA_real_
  out <- embeddings * s
  dim(out) <- dim(embeddings)
  out
}

# Per-zone probe of projected points, as in the Python accessor's probe().
probe_zone <- function(gt, zone, e, n, year, search_px, array, n_bands) {
  k <- length(e)
  status <- rep(OUTSIDE, k)
  values <- matrix(NA_real_, k, n_bands)
  ok <- is.finite(e) & is.finite(n) & zone$width > 0 & zone$height > 0
  xi <- yi <- rep(NA_real_, k)
  xi[ok] <- nearest_index(e[ok], zone$ox, zone$px, zone$width)
  yi[ok] <- nearest_index(n[ok], zone$oy, -zone$px, zone$height)
  inside <- ok & abs(zone_x(zone, xi) - e) <= zone$px &
    abs(zone_y(zone, yi) - n) <= zone$px
  idx <- which(inside)
  if (!length(idx)) return(list(status = status, values = values))

  ti <- zone_time_index(zone, year)
  scales_reader <- zone_reader(gt, zone, "scales")
  emb_reader <- zone_reader(gt, zone, array)
  rows <- cbind(yi[idx], xi[idx])
  load_indices(list(
    list(reader = scales_reader,
         keys = shard_location(scales_reader, floor(cbind(ti, rows) /
           rep(scales_reader$inner, each = length(idx))))$keys),
    list(reader = emb_reader,
         keys = shard_location(emb_reader, floor(cbind(ti, 0, rows) /
           rep(emb_reader$inner, each = length(idx))))$keys)
  ))
  centre <- read_points(scales_reader, cbind(ti, yi[idx], xi[idx]))[, 1L]
  water <- is.nan(centre) & !zone$source_nodata
  valid <- is.finite(centre)
  status[idx[water]] <- WATER
  bx <- xi
  by <- yi
  scale <- rep(NA_real_, k)
  scale[idx[valid]] <- centre[valid]

  r <- search_px
  for (p in idx[!valid & !water]) {
    x0 <- max(0, xi[p] - r)
    x1 <- min(zone$width, xi[p] + r + 1)
    y0 <- max(0, yi[p] - r)
    y1 <- min(zone$height, yi[p] + r + 1)
    win <- read_block(scales_reader, c(ti, y0, x0), c(ti + 1, y1, x1))
    win <- matrix(win, nrow = x1 - x0)
    cand <- which(t(is.finite(win)), arr.ind = TRUE)
    if (!nrow(cand)) {
      status[p] <- NODATA
      next
    }
    cand <- cand[order(cand[, 1L], cand[, 2L]), , drop = FALSE]
    d <- (cand[, 1L] - 1 - (yi[p] - y0))^2 + (cand[, 2L] - 1 - (xi[p] - x0))^2
    best <- which.min(d)
    by[p] <- y0 + cand[best, 1L] - 1
    bx[p] <- x0 + cand[best, 2L] - 1
    scale[p] <- win[cand[best, 2L], cand[best, 1L]]
  }

  good <- which(is.finite(scale))
  if (length(good)) {
    emb <- read_points(emb_reader, cbind(ti, 0, by[good], bx[good]), whole = 2L)
    values[good, ] <- emb * scale[good]
    status[good] <- VALID
  }
  list(status = status, values = values)
}
