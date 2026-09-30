test_that("seam neighbours cover half a tile each side", {
  nb <- function(lon) {
    x <- seam_neighbours(lon)
    as.integer(stats::na.omit(x[1, ]))
  }
  expect_identical(nb(3.0), integer(0))
  expect_identical(nb(138.03), 53L)
  expect_identical(nb(137.97), 54L)
  expect_identical(nb(138.2), integer(0))
  expect_identical(nb(138.0), 53L)
  expect_identical(nb(-179.97), 60L)
  expect_identical(nb(179.97), 1L)
})

test_that("zone routing follows the longitude rule", {
  expect_identical(zone_for_lon(c(2.35, -120.5, -180, 180)), c(31L, 10L, 1L, 60L))
})

test_that("zones spanned walk the short way round", {
  expect_identical(zones_spanned(c(-0.1, 0.1), 0), c(30L, 31L))
  expect_identical(zones_spanned(c(179.9, -179.9), 180), c(60L, 1L))
  expect_identical(zones_spanned(c(1, 2), 1.5), 31L)
})

test_that("store URLs resolve published versions", {
  base <- "https://data.source.coop/tessera/tessera/zarr/"
  expect_identical(tessera_store_url(), paste0(base, "v1.1-dclimate"))
  expect_identical(tessera_store_url("1.0"), paste0(base, "v1"))
  expect_identical(tessera_store_url("v1"), paste0(base, "v1"))
  expect_identical(tessera_store_url("1.1"), paste0(base, "v1.1-dclimate"))
  expect_identical(tessera_store_url("v1.1", "cambridge"), paste0(base, "v1.1"))
  expect_identical(tessera_store_url("v2"), paste0(base, "v2-2B-L~beta1"))
  expect_identical(tessera_store_url("v2", "2B-L~beta2"), paste0(base, "v2-2B-L~beta2"))
  expect_identical(tessera_store_url("v1.1-dclimate"), paste0(base, "v1.1-dclimate"))
  expect_identical(tessera_store_url("v3"), paste0(base, "v3"))
  expect_identical(tessera_store_url("v3", "x"), paste0(base, "v3-x"))
  expect_error(tessera_store_url("v1.1", "nope"), "available: dclimate, cambridge")
  expect_error(tessera_store_url("v1.1-dclimate", "cambridge"), "needs a version")
})

test_that("published datasets resolve from their store URLs", {
  ds <- tessera_datasets()
  expect_identical(ds$default, c(TRUE, TRUE, FALSE, TRUE, FALSE))
  for (i in seq_len(nrow(ds))) {
    expect_identical(dataset_for_location(paste0(ds$url[i], "/")),
                     list(version = ds$version[i], variant = ds$variant[i]))
  }
  expect_null(dataset_for_location("/data/tessera.zarr"))
  expect_identical(dataset_version_path(c("1.0", "1.1", "2.0")), c("v1", "v1.1", "v2"))
  expect_identical(dataset_tags(list(version = "1.1", variant = "dclimate")),
                   c(TESSERA_DATASET_VERSION = "1.1",
                     TESSERA_DATASET_VERSION_PATH = "v1.1",
                     TESSERA_DATASET_VARIANT = "dclimate"))
})

test_that("incomplete years warn only for the hemisphere read", {
  groups <- list("31N" = list(years_complete = list(2023L, 2024L)),
                 "31S" = list(years_complete = list(2024L)),
                 other = list(crs = "EPSG:32631"))
  inc <- incomplete_years(groups, c(2023L, 2024L))
  expect_identical(inc, list("31S" = 2023L))
  zone <- list(oy = 1000, px = 10, incomplete = inc)
  expect_silent(warn_incomplete(zone, 2023, 0, 50))
  expect_warning(warn_incomplete(zone, 2023, 50, 150), "31S has no complete 2023",
                 class = "geotessera_warning")
  expect_silent(warn_incomplete(zone, 2024, 50, 150))
})

test_that("cache keys separate store locations", {
  expect_identical(store_cache_key(tessera_store_url("v2")), "v2-2B-L_beta1")
  a <- store_cache_key("https://example.org/a/v1")
  b <- store_cache_key("https://example.com/a/v1")
  expect_false(identical(a, b))
  expect_match(a, "^example.org_a_v1-[0-9a-f]{12}$")
})

test_that("range headers use fixed notation", {
  expect_identical(range_header(c(100000, 300000)), "bytes=100000-299999")
  expect_identical(range_header(-262148), "bytes=-262148")
  expect_identical(range_header(5e9), "bytes=5000000000-")
  body <- as.raw(0:9)
  expect_identical(slice_range(body, c(2, 5)), as.raw(2:4))
  expect_identical(slice_range(body, -3), as.raw(7:9))
})

test_that("nearest pixels break ties towards the lower index", {
  expect_identical(nearest_index(c(10, 14.9, 15, 15.1), 0, 10, 5), c(0, 1, 1, 1))
  expect_identical(nearest_index(c(-100, 100), 0, 10, 5), c(0, 4))
  expect_identical(nearest_index(c(95, 90, 85), 100, -10, 5), c(0, 0, 1))
})

test_that("points accept matrices, vectors and data frames", {
  m <- cbind(c(1, 2), c(3, 4))
  expect_identical(as_lonlat(m), m)
  expect_identical(as_lonlat(c(1, 3)), matrix(c(1, 3), 1))
  expect_identical(as_lonlat(data.frame(lat = c(3, 4), lon = c(1, 2))), m)
  expect_error(as_lonlat("a"), "longitude and latitude")
})

test_that("dequantizing masks non-finite scales", {
  emb <- array(1:8, c(2, 2, 2))
  scales <- matrix(c(0.5, NaN, Inf, 2), 2)
  out <- dequantize_array(emb, scales)
  expect_equal(out[, , 1], matrix(c(0.5, NA, NA, 8), 2))
  expect_equal(out[, , 2], matrix(c(2.5, NA, NA, 16), 2))
})

test_that("data frames without lon and lat columns are ambiguous", {
  expect_identical(as_lonlat(data.frame(a = c(1, 2), b = c(3, 4))),
                   cbind(c(1, 2), c(3, 4)))
  expect_error(as_lonlat(data.frame(id = 1:2, a = c(1, 2), b = c(3, 4))),
               "`lon` and `lat` columns", class = "geotessera_error")
})

test_that("bounds accept vectors, extents and spatial objects", {
  bbox <- c(0.004, 52.005, 0.013, 52.009)
  lonlat <- function(b) as_bounds(b)$lonlat
  expect_identical(lonlat(bbox), bbox)
  expect_identical(lonlat(c(ymax = 52.009, xmin = 0.004, xmax = 0.013,
                            ymin = 52.005)), bbox)
  expect_identical(lonlat(terra::ext(0.004, 0.013, 52.005, 52.009)), bbox)
  v <- terra::vect(rbind(c(0.004, 52.005), c(0.013, 52.009)), crs = "EPSG:4326")
  expect_equal(lonlat(v), bbox)

  # Projected bounds keep their CRS and box; their WGS84 bounds enclose them.
  utm <- terra::project(terra::as.polygons(terra::ext(v), crs = "EPSG:4326"),
                        "EPSG:32631")
  b <- as_bounds(utm)
  expect_identical(b$box, unname(as.vector(terra::ext(utm))[c(1, 3, 2, 4)]))
  expect_identical(bounds_envelope(b, "EPSG:32631"), b$box)
  expect_true(all(b$lonlat[1:2] <= bbox[1:2] + 1e-9) &&
                all(b$lonlat[3:4] >= bbox[3:4] - 1e-9))

  expect_error(as_bounds(terra::vect(rbind(c(0, 0), c(1, 1)))), "declare its CRS")
  expect_error(as_bounds(c(179, 0, -179, 1)), "antimeridian", class = "geotessera_error")
  expect_error(as_bounds(1:3), "four finite coordinates")
})

test_that("sf bounds keep their CRS", {
  skip_if_not_installed("sf")
  bbox <- c(0.004, 52.005, 0.013, 52.009)
  b <- sf::st_bbox(c(xmin = 0.004, ymin = 52.005, xmax = 0.013, ymax = 52.009),
                   crs = sf::st_crs(4326))
  expect_equal(as_bounds(b)$lonlat, bbox)
  expect_equal(as_bounds(sf::st_as_sfc(b))$lonlat, bbox)
  utm <- sf::st_transform(sf::st_as_sfc(b), 32631)
  expect_equal(as_bounds(utm)$box, unname(as.numeric(sf::st_bbox(utm))))
})

test_that("retry backoff leaves the random number stream alone", {
  set.seed(1)
  before <- .Random.seed
  j <- retry_jitter()
  expect_identical(.Random.seed, before)
  expect_true(j >= 0.5 && j < 1)
})

test_that("disk caches persist between sessions", {
  dir <- withr::local_tempdir()
  a <- new_cache(1024^2, dir)
  cache_put(a, "index:utm31/embeddings/c/0/0/0/0", list(offset = 1, length = 2))
  cache_put(a, "utm31/embeddings/c/0/0/0/0#3", as.raw(1:4))
  b <- new_cache(1024^2, dir)
  expect_identical(cache_get(b, "utm31/embeddings/c/0/0/0/0#3"), as.raw(1:4))
  expect_identical(cache_get(b, "index:utm31/embeddings/c/0/0/0/0"),
                   list(offset = 1, length = 2))
  expect_null(cache_get(b, "absent"))
  expect_null(cache_get(new_cache(1024^2), "utm31/embeddings/c/0/0/0/0#3"))
})
