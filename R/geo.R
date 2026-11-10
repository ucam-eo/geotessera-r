# Zone routing and coordinate helpers.

wgs84 <- "EPSG:4326"

# Tiles are 0.1 degrees and belong to the zone holding their centre, so a
# point within half a tile of a seam may be covered by the neighbouring
# zone. This doubles that bound to absorb curved tile footprints.
seam_degrees <- 0.1

zone_for_lon <- function(lon) {
  as.integer(pmax(1, pmin(60, floor((lon + 180) / 6) + 1)))
}

# Neighbouring zones to try after a point's own zone, as columns `west`
# and `east`; NA away from a seam.
seam_neighbours <- function(lon) {
  z <- zone_for_lon(lon)
  frac <- (lon + 180) %% 6
  west <- ifelse(frac <= seam_degrees, ifelse(z == 1L, 60L, z - 1L), NA_integer_)
  east <- ifelse(frac >= 6 - seam_degrees, ifelse(z == 60L, 1L, z + 1L), NA_integer_)
  cbind(west = west, east = east)
}

# The contiguous run of zones covering `lons`, walked the short way round
# so that a run across the antimeridian is c(60, 1).
zones_spanned <- function(lons, centre_lon) {
  zc <- zone_for_lon(centre_lon)
  d <- (zone_for_lon(lons) - zc) %% 60
  offsets <- ifelse(d > 30, d - 60, d)
  (zc - 1L + seq(min(offsets), max(offsets))) %% 60L + 1L
}

same_crs <- function(a, b) {
  identical(a, b) || isTRUE(terra::same.crs(a, b))
}

project_xy <- function(x, y, from, to) {
  if (!length(x) || same_crs(from, to)) return(list(x = x, y = y))
  xy <- terra::project(cbind(as.numeric(x), as.numeric(y)), from = from, to = to)
  list(x = xy[, 1L], y = xy[, 2L])
}

# Points along the edges of a box c(xmin, ymin, xmax, ymax), so that its
# curved outline survives projection.
box_outline <- function(box, densify = 21L) {
  t <- seq(0, 1, length.out = densify + 2L)
  w <- box[1L]; s <- box[2L]; e <- box[3L]; n <- box[4L]
  cbind(c(w + (e - w) * t, w + (e - w) * t, rep(w, length(t)), rep(e, length(t))),
        c(rep(s, length(t)), rep(n, length(t)), s + (n - s) * t, s + (n - s) * t))
}

# The extent in `to` of a box in `from`, including its curved edges. A box
# already in `to` is returned unchanged.
box_envelope <- function(box, from, to) {
  if (same_crs(from, to)) return(box)
  outline <- box_outline(box)
  xy <- project_xy(outline[, 1L], outline[, 2L], from, to)
  c(min(xy$x), min(xy$y), max(xy$x), max(xy$y))
}

# The projected extent of a lon/lat bbox.
utm_envelope <- function(bbox, crs) box_envelope(bbox, wgs84, crs)

# Bounds from any accepted form, as `box` c(xmin, ymin, xmax, ymax) in its
# own `crs`, and `lonlat`, the WGS84 bounds c(west, south, east, north)
# used to choose zones. Windows are computed from `box` with at most one
# transform, so bounds given in a zone's CRS select exactly the pixels
# whose centres they enclose.
as_bounds <- function(bbox) {
  corners <- c("xmin", "ymin", "xmax", "ymax")
  crs <- wgs84
  if (inherits(bbox, c("sf", "sfc", "bbox"))) {
    check_installed("sf", "sf bounds")
    if (is.na(sf::st_crs(bbox))) abort("`bbox` must declare its CRS")
    crs <- sf::st_crs(bbox)$wkt
    bbox <- unclass(sf::st_bbox(bbox))
  } else if (inherits(bbox, c("SpatVector", "SpatRaster"))) {
    crs <- terra::crs(bbox)
    if (!nzchar(crs)) abort("`bbox` must declare its CRS")
    bbox <- as.vector(terra::ext(bbox))
  } else if (inherits(bbox, "SpatExtent")) {
    bbox <- as.vector(bbox)
  }
  if (all(corners %in% names(bbox))) bbox <- bbox[corners]
  bbox <- unname(bbox)
  if (!is.numeric(bbox) || length(bbox) != 4L || !all(is.finite(bbox))) {
    abort("`bbox` must be four finite coordinates c(west, south, east, north)")
  }
  box <- as.numeric(bbox)
  lonlat <- box_envelope(box, crs, wgs84)
  if (!(box[1L] < box[3L] && box[2L] < box[4L] &&
        -180 <= lonlat[1L] && lonlat[3L] <= 180 &&
        -90 <= lonlat[2L] && lonlat[4L] <= 90)) {
    abort("`bbox` must have positive area within WGS84 bounds; ",
          "split bounds across the antimeridian into two requests")
  }
  list(box = box, crs = crs, lonlat = lonlat)
}

# The extent of bounds in `crs`.
bounds_envelope <- function(bounds, crs) box_envelope(bounds$box, bounds$crs, crs)

# Points from any accepted form, as `xy`, coordinates in their own `crs`,
# and `lonlat`, their WGS84 coordinates used to choose zones.
as_points <- function(points) {
  if (inherits(points, c("sf", "sfc"))) {
    check_installed("sf", "sf points")
    if (is.na(sf::st_crs(points))) {
      abort("Points must declare their CRS")
    }
    if (!all(sf::st_geometry_type(points) == "POINT")) {
      abort("Points must be POINT geometries")
    }
    xy <- unname(sf::st_coordinates(points)[, 1:2, drop = FALSE])
    crs <- sf::st_crs(points)$wkt
  } else if (inherits(points, "SpatVector")) {
    crs <- terra::crs(points)
    if (!nzchar(crs)) {
      abort("Points must declare their CRS")
    }
    xy <- unname(terra::crds(points))
    if (terra::geomtype(points) != "points" || nrow(xy) != nrow(points)) {
      abort("Points must be single point geometries")
    }
  } else {
    xy <- as_lonlat(points)
    crs <- wgs84
  }
  ll <- project_xy(xy[, 1L], xy[, 2L], crs, wgs84)
  list(xy = xy, crs = crs, lonlat = cbind(ll$x, ll$y))
}

# Longitude and latitude columns from points given as coordinates.
as_lonlat <- function(points) {
  if (is.data.frame(points)) {
    cols <- tolower(names(points))
    lon <- match(c("lon", "longitude", "x"), cols)
    lat <- match(c("lat", "latitude", "y"), cols)
    lon <- lon[!is.na(lon)][1L]
    lat <- lat[!is.na(lat)][1L]
    if (is.na(lon) || is.na(lat)) {
      if (ncol(points) != 2L) {
        abort("A data frame of points needs `lon` and `lat` columns")
      }
      lon <- 1L
      lat <- 2L
    }
    points <- cbind(as.numeric(points[[lon]]), as.numeric(points[[lat]]))
  }
  if (is.numeric(points) && is.null(dim(points)) && length(points) == 2L) {
    points <- matrix(points, ncol = 2L)
  }
  if (!is.numeric(points) || length(dim(points)) != 2L || ncol(points) != 2L) {
    abort("Points must be longitude and latitude pairs")
  }
  if (any(abs(points[, 1L]) > 180 | abs(points[, 2L]) > 90, na.rm = TRUE)) {
    abort("Points must be WGS84 longitudes and latitudes; give projected ",
          "coordinates as an sf or SpatVector layer with their CRS")
  }
  unname(points)
}
