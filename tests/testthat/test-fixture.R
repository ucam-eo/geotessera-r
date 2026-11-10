# Outputs of the Python reader on the same store; see data-raw/make-fixture.py.

test_that("the store opens with its metadata", {
  gt <- fixture_store()
  expect_s3_class(gt, "tessera")
  expect_identical(gt$years, c(2023L, 2024L))
  expect_identical(gt$n_bands, 4L)
  expect_identical(gt$zone_ids, c(30L, 31L))
  expect_identical(unname(gt$depths[c("2", "4")]), c("embeddings_d2", "embeddings"))
  expect_output(print(gt), "Years  : 2023, 2024")
})

test_that("local stores are recorded by absolute path", {
  withr::local_dir(dirname(tessera_example()))
  gt <- tessera_open("example.zarr", cache_dir = "cache")
  expect_identical(gt$url, normalizePath(tessera_example()))
  expect_identical(gt$cache_dir, normalizePath("cache", mustWork = FALSE))
})

test_that("a store without consolidated metadata or zones is an error", {
  root <- jsonlite::read_json(file.path(tessera_example(), "zarr.json"))
  write_store <- function(meta) {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    jsonlite::write_json(meta, file.path(dir, "zarr.json"), auto_unbox = TRUE,
                         null = "null")
    dir
  }
  bare <- root
  bare$consolidated_metadata <- NULL
  expect_error(tessera_open(write_store(bare)), "has no consolidated metadata")
  zoneless <- root
  keys <- names(zoneless$consolidated_metadata$metadata)
  zoneless$consolidated_metadata$metadata <-
    zoneless$consolidated_metadata$metadata[!grepl("^utm", keys)]
  zoneless$consolidated_metadata$metadata$other <- list(node_type = "group",
                                                        zarr_format = 3L)
  expect_error(tessera_open(write_store(zoneless)), "No UTM zone groups")
})

test_that("probe and sample match the Python reader", {
  gt <- fixture_store()
  g <- golden()
  pts <- golden_points(g)
  for (year in c(2023, 2024)) {
    res <- tessera_probe(gt, pts, year)
    expect_s3_class(res$status, "factor")
    expect_identical(as.character(res$status), g[[paste0("probe_", year)]])
    expect_values(band_matrix(res), g[[paste0("sample_", year)]])
    expect_identical(names(res), c("lon", "lat", paste0("band_", 1:4), "status"))
    plain <- tessera_sample(gt, pts, year, search_px = 0, cross_zone = FALSE)
    expect_values(band_matrix(plain), g[[paste0("sample_nosearch_", year)]])
  }
  expect_values(band_matrix(tessera_sample(gt, pts, 2024, depth = 2)), g$sample_d2)
})

test_that("points outside coverage or non-finite are NA", {
  gt <- fixture_store()
  res <- tessera_probe(gt, rbind(c(NA, 52), c(100, 10)), 2024)
  expect_identical(as.character(res$status), c("outside", "outside"))
  expect_true(all(is.na(band_matrix(res))))
  expect_equal(dim(tessera_sample(gt, matrix(numeric(0), 0, 2), 2024)), c(0L, 6L))
  expect_error(tessera_sample(gt, data.frame(x = 703500, y = 5766500), 2024),
               "Points must be WGS84 longitudes and latitudes")
  expect_error(tessera_sample(gt, c(0, 91), 2024), "WGS84")
})

test_that("invalid years and depths are reported", {
  gt <- fixture_store()
  pt <- golden_points(golden())[1, ]
  expect_error(tessera_sample(gt, pt, 2019), "Year 2019 is unavailable; available: 2023, 2024")
  expect_error(tessera_sample(gt, c(10, 10), 2019), "Year 2019 is unavailable")
  expect_error(tessera_sample(gt, pt, 2024.5), "`year` must be a single whole number")
  expect_error(tessera_sample(gt, pt, c(2023, 2024)), "`year` must be a single whole number")
  bbox <- c(0.004, 52.005, 0.013, 52.009)
  expect_error(tessera_region(gt, bbox, 2019), "Year 2019 is unavailable")
  expect_error(tessera_patch(gt, pt, 2019, size_px = 4), "Year 2019 is unavailable")
  expect_error(tessera_export_geotiffs(gt, bbox, 2019, tempdir(), dry_run = TRUE),
               "Year 2019 is unavailable")
  expect_error(tessera_sample(gt, pt, 2024, depth = 3), "depth 3 not in this store; available: 2, 4")
})

test_that("regions match the Python reader", {
  gt <- fixture_store()
  g <- golden()$region
  r <- tessera_region(gt, g$bbox, 2024)
  expect_s4_class(r, "SpatRaster")
  expect_equal(dim(r), as.integer(g$shape))
  tr <- g$transform
  expect_equal(as.vector(terra::ext(r)),
               c(xmin = tr[3], xmax = tr[3] + g$shape[2] * tr[1],
                 ymin = tr[6] + g$shape[1] * tr[5], ymax = tr[6]))
  expect_identical(terra::crs(r, describe = TRUE)$code, "32631")
  expect_values(terra::as.array(r), g$values)

  q <- tessera_region_quantized(gt, g$bbox, 2024)
  expect_s4_class(q, "SpatRaster")
  expect_identical(names(q), c(paste0("band_", 1:4), "scale"))
  expect_equal(unname(terra::as.array(q[[1:4]])), g$quantized, ignore_attr = TRUE)
  expect_values(terra::as.array(q[["scale"]])[, , 1], g$scales)
  expect_values(terra::as.array(tessera_dequantize(q)), g$values)
  expect_equal(as.vector(terra::ext(q)), as.vector(terra::ext(r)))
  expect_error(tessera_dequantize(r), "tessera_region_quantized")

  sub <- tessera_region(gt, g$bbox, 2024, bands = c(3, 1))
  expect_identical(names(sub), c("band_3", "band_1"))
  expect_values(terra::as.array(sub), g$values[, , c(3, 1)])
})

test_that("region strips cover the region top to bottom", {
  gt <- fixture_store()
  g <- golden()$region
  strips <- tessera_map_region(gt, g$bbox, 2024, function(s) s, strip_rows = 16)
  expect_length(strips, length(g$strips))
  for (i in seq_along(strips)) {
    expect_equal(dim(strips[[i]]), as.integer(g$strips[[i]][[1]]))
    expect_equal(terra::ymax(strips[[i]]), g$strips[[i]][[2]][6])
  }
  merged <- do.call(terra::merge, unname(strips))
  expect_values(terra::as.array(merged), g$values)
})

test_that("patches match the Python reader", {
  gt <- fixture_store()
  g <- golden()
  for (key in c("patch_native", "patch_edge", "patch_seam")) {
    p <- g[[key]]
    # Python merges seam patches onto a transverse Mercator grid centred on
    # the point; the same grid here exercises the resampling path.
    dst <- if (key == "patch_seam") sprintf(
      "+proj=tmerc +lat_0=0 +lon_0=%.8f +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs",
      p$lon)
    r <- tessera_patch(gt, c(p$lon, p$lat), 2023, p$size, dst_crs = dst)
    expect_equal(dim(r), c(p$size, p$size, 4L))
    expect_equal(terra::xmin(r), p$transform[3], tolerance = 1e-6)
    expect_equal(terra::ymax(r), p$transform[6], tolerance = 1e-6)
    expect_values(terra::as.array(r), p$values)
    if (p$native) expect_identical(terra::crs(r, describe = TRUE)$code, "32631")
  }
})

test_that("patches honour a destination CRS", {
  gt <- fixture_store()
  p <- golden()$patch_native
  r <- tessera_patch(gt, c(p$lon, p$lat), 2023, 8, dst_crs = "EPSG:3857")
  expect_identical(terra::crs(r, describe = TRUE)$code, "3857")
  expect_error(tessera_patch(gt, c(p$lon, p$lat), 2023, 8, dst_crs = "EPSG:4326"),
               "projected CRS with metre units")
  expect_error(tessera_patch(gt, c(p$lon, p$lat), 2023, 0), "positive integer")
})

test_that("export estimates match the Python exporter", {
  gt <- fixture_store()
  g <- golden()$export
  est <- tessera_export_geotiffs(gt, g$bbox, 2024, withr::local_tempfile(), dry_run = TRUE)
  expect_equal(est, g$estimates, ignore_attr = TRUE)
})

test_that("exported GeoTIFFs hold the native pixels and metadata", {
  gt <- fixture_store()
  g <- golden()$export
  dir <- withr::local_tempfile()
  calls <- 0
  files <- tessera_export_geotiffs(gt, g$bbox, 2024, dir, bands = c(2, 4),
                                   strip_rows = 7, progress = function(...) calls <<- calls + 1)
  expect_identical(basename(files), c("tessera_2024_utm30.tif", "tessera_2024_utm31.tif"))
  expect_identical(sort(list.files(dir)), basename(files))
  expect_identical(calls, 8)

  r <- terra::rast(files[2])
  expect_equal(dim(r), c(g$estimates$height[2], g$estimates$width[2], 2))
  expect_identical(names(r), c("band_2", "band_4"))
  tags <- terra::metags(r)
  tags <- stats::setNames(tags$value, tags$name)
  expect_identical(unname(tags[c("TESSERA_YEAR", "TESSERA_BANDS", "TESSERA_MODEL")]),
                   c("2024", "2,4", "https://geotessera.org/model/test"))

  zone <- zone_open(gt, 31)
  e <- terra::ext(r)
  win <- zone_window(zone, c(e$xmin + 5, e$ymin + 5, e$xmax - 5, e$ymax - 5))
  q <- read_window(gt, zone, "embeddings", 1L, win$x0, win$x1, win$y0, win$y1, c(2L, 4L))
  expected <- aperm(dequantize_array(q$embeddings, q$scales), c(2, 1, 3))
  expect_values(terra::as.array(r), expected)
})

test_that("SpatVector points are projected to WGS84", {
  gt <- fixture_store()
  g <- golden()
  pts <- golden_points(g)[1:7, ]
  v <- terra::project(terra::vect(pts, crs = "EPSG:4326"), "EPSG:32631")
  out <- tessera_sample(gt, v, 2024)
  expect_s4_class(out, "SpatVector")
  expect_values(band_matrix(terra::values(out)), g$sample_2024[1:7, ])
})

test_that("regions accept spatial bounds in any CRS", {
  gt <- fixture_store()
  g <- golden()$region
  b <- g$bbox
  e <- terra::ext(b[1], b[3], b[2], b[4])
  expect_values(terra::as.array(tessera_region(gt, e, 2024)), g$values)
  expect_error(tessera_region(gt, c(b[3], b[2], b[1], b[4]), 2024),
               "positive area", class = "geotessera_error")
})

test_that("sf points are projected to WGS84", {
  skip_if_not_installed("sf")
  gt <- fixture_store()
  g <- golden()
  pts <- golden_points(g)[1:7, ]
  sfp <- sf::st_transform(sf::st_as_sf(as.data.frame(pts), coords = 1:2, crs = 4326),
                          32631)
  out <- tessera_sample(gt, sfp, 2024)
  expect_s3_class(out, "sf")
  expect_values(band_matrix(sf::st_drop_geometry(out)), g$sample_2024[1:7, ])
  expect_s3_class(tessera_sample(gt, sf::st_geometry(sfp), 2024), "sf")
})

test_that("errors are classed", {
  gt <- fixture_store()
  expect_error(tessera_sample(gt, c(0, 52), 2019), class = "geotessera_error")
  expect_error(tessera_region("gt", c(0, 52, 1, 53), 2024), class = "geotessera_error")
})

test_that("exports return their paths invisibly", {
  gt <- fixture_store()
  g <- golden()$export
  dir <- withr::local_tempfile()
  expect_invisible(tessera_export_geotiffs(gt, g$bbox, 2024, dir, bands = 1))
  expect_error(tessera_export_geotiffs(gt, g$bbox, 2024, dir, compress = "jpeg"),
               "should be one of")
})

test_that("samples keep the class and columns of their points", {
  gt <- fixture_store()
  g <- golden()
  pts <- golden_points(g)[1:3, ]
  df <- data.frame(id = c("a", "b", "c"), longitude = pts[, 1], latitude = pts[, 2])
  out <- tessera_sample(gt, df, 2024)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c(names(df), paste0("band_", 1:4)))
  expect_identical(out$id, df$id)
  expect_values(band_matrix(out), g$sample_2024[1:3, ])

  one <- tessera_sample(gt, pts[1, ], 2024)
  expect_identical(names(one), c("lon", "lat", paste0("band_", 1:4)))

  df$band_1 <- 0
  expect_error(tessera_sample(gt, df, 2024), "already has columns named `band_1`")
})

test_that("seam patches keep the native grid of the point's zone", {
  gt <- fixture_store()
  p <- golden()$patch_seam
  r <- tessera_patch(gt, c(p$lon, p$lat), 2023, p$size)
  zone <- zone_open(gt, 31)
  expect_true(same_crs(terra::crs(r), zone$crs))
  off <- grid_offsets(zone, terra::ext(r), zone$crs)
  expect_false(is.null(off))

  a <- terra::as.array(r)
  native <- terra::as.array(patch_native(gt, zone, off$x0, off$y0, 2023, p$size,
                                         list(name = "embeddings"), NULL, 1:4))
  cells <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  ll <- terra::project(cells, from = zone$crs, to = "EPSG:4326")
  owner <- matrix(zone_for_lon(ll[, 1]), p$size, p$size, byrow = TRUE)
  expect_true(any(owner == 30) && any(owner == 31))

  # Pixels of zone 31 are its stored pixels, unresampled.
  mine <- owner == 31 & !is.na(native[, , 1])
  expect_identical(a[, , 1][mine], native[, , 1][mine])

  # Pixels of zone 30 are its nearest pixels.
  theirs <- which(t(owner == 30 & !is.na(a[, , 1])))
  s <- tessera_sample(gt, ll[theirs, , drop = FALSE], 2023, search_px = 0,
                      cross_zone = FALSE)
  got <- terra::values(r)[theirs, ]
  # GDAL warps in float32.
  expect_equal(unname(got), unname(band_matrix(s)), tolerance = 1e-6)

  # dst_crs naming the zone's CRS is the same native grid.
  same <- tessera_patch(gt, c(p$lon, p$lat), 2023, p$size, dst_crs = "EPSG:32631")
  expect_identical(terra::values(same), terra::values(r))
})

test_that("bounds in a zone's CRS select exactly the pixels they enclose", {
  gt <- fixture_store()
  zone <- zone_open(gt, 31)
  ex <- terra::ext(zone$ox + 200, zone$ox + 600, zone$oy - 900, zone$oy - 600)
  tmpl <- terra::rast(ex, res = 10, crs = zone$crs)
  r <- tessera_region(gt, tmpl, 2024)
  expect_equal(dim(r)[1:2], dim(tmpl)[1:2])
  expect_equal(as.vector(terra::ext(r)), as.vector(ex))
  expect_true(same_crs(terra::crs(r), zone$crs))

  # Across the seam, zone 31's grid is still chosen for zone 31 bounds.
  seam <- terra::ext(zone$ox - 300, zone$ox + 300, zone$oy - 900, zone$oy - 600)
  r2 <- suppressWarnings(tessera_region(gt, terra::rast(seam, res = 10, crs = zone$crs),
                                        2024))
  expect_true(same_crs(terra::crs(r2), zone$crs))

  est <- tessera_export_geotiffs(gt, tmpl, 2024, tempdir(), dry_run = TRUE)
  expect_identical(est$zone, 31L)
  expect_equal(c(est$width, est$height), c(40, 30))
})

test_that("points in a zone's CRS are sampled without a round trip", {
  gt <- fixture_store()
  zone <- zone_open(gt, 31)
  # Exactly on a pixel edge; ties choose the lower index.
  e <- zone$ox + 20 * 10
  n <- zone$oy - 60 * 10 - 5
  v <- terra::vect(cbind(e, n), crs = zone$crs)
  q <- read_window(gt, zone, "embeddings", 1L, 19, 20, 60, 61)
  expected <- dequantize_array(q$embeddings, q$scales)[1, 1, ]
  got <- band_matrix(terra::values(tessera_sample(gt, v, 2024, search_px = 0)))
  expect_equal(unname(got[1, ]), expected)
})

test_that("exports record their dataset and refuse to mix datasets", {
  gt <- fixture_store()
  expect_null(gt$dataset)
  g <- golden()$export
  dir <- withr::local_tempfile()
  tessera_export_geotiffs(gt, g$bbox, 2024, dir, bands = 1)
  expect_false(file.exists(file.path(dir, "tessera_metadata.json")))

  gt$dataset <- list(version = "1.1", variant = "dclimate")
  expect_output(print(gt), "Dataset: v1.1 dclimate")
  files <- tessera_export_geotiffs(gt, g$bbox, 2024, dir, bands = 1)
  meta <- jsonlite::read_json(file.path(dir, "tessera_metadata.json"))
  expect_identical(meta[c("dataset_version", "dataset_version_path", "dataset_variant",
                          "format", "year")],
                   list(dataset_version = "1.1", dataset_version_path = "v1.1",
                        dataset_variant = "dclimate", format = "tiff", year = 2024L))
  tags <- terra::metags(terra::rast(files[1]))
  tags <- stats::setNames(tags$value, tags$name)
  expect_identical(unname(tags[c("TESSERA_DATASET_VERSION", "TESSERA_DATASET_VARIANT")]),
                   c("1.1", "dclimate"))
  # The same dataset may export again into the directory.
  expect_no_error(tessera_export_geotiffs(gt, g$bbox, 2023, dir, bands = 1))

  gt$dataset <- list(version = "1.1", variant = "cambridge")
  expect_error(tessera_export_geotiffs(gt, g$bbox, 2024, dir, bands = 1),
               "holds v1.1 dclimate embeddings", class = "geotessera_dataset_error")
  # A dry run writes nothing and so is not checked.
  expect_s3_class(tessera_export_geotiffs(gt, g$bbox, 2024, dir, dry_run = TRUE),
                  "data.frame")
})

test_that("samples select bands and validate their options", {
  gt <- fixture_store()
  g <- golden()
  pts <- golden_points(g)
  out <- tessera_sample(gt, pts, 2024, bands = c(3, 1))
  expect_identical(names(out), c("lon", "lat", "band_3", "band_1"))
  expect_values(band_matrix(out), g$sample_2024[, c(3, 1)])
  expect_error(tessera_sample(gt, pts, 2024, bands = 5), "no greater than 4")
  expect_error(tessera_sample(gt, pts, 2024, search_px = -1), "non-negative integer")
  expect_error(tessera_sample(gt, pts, 2024, cross_zone = NA), "TRUE or FALSE")
})

test_that("every region reader warns of bounds across zones", {
  gt <- fixture_store()
  bbox <- golden()$export$bbox
  expect_warning(tessera_region(gt, bbox, 2024), "spans UTM zones 30, 31")
  expect_warning(tessera_region_quantized(gt, bbox, 2024), "spans UTM zones",
                 class = "geotessera_warning")
  expect_warning(tessera_map_region(gt, bbox, 2024, dim), "spans UTM zones")
})

test_that("strips stay within their memory bound", {
  expect_identical(strip_height(100, 4, 16L), 16L)
  expect_identical(strip_height(2000, 128, 512L), 13L)
  expect_identical(strip_height(1e9, 128, 512L), 1L)
})

test_that("stores and output directories are validated", {
  expect_error(tessera_open(NA_character_), "URL or a path")
  expect_error(tessera_open("s3://bucket/store"), "tessera_datasets",
               class = "geotessera_error")
  gt <- fixture_store()
  expect_error(tessera_export_geotiffs(gt, golden()$export$bbox, 2024, c("a", "b")),
               "single path")
})

test_that("patch centres accept any single point", {
  gt <- fixture_store()
  p <- golden()$patch_native
  ref <- tessera_patch(gt, c(p$lon, p$lat), 2023, p$size)
  forms <- list(
    data.frame(lon = p$lon, lat = p$lat),
    cbind(p$lon, p$lat),
    terra::vect(cbind(p$lon, p$lat), crs = "EPSG:4326")
  )
  for (point in forms) {
    expect_identical(terra::values(tessera_patch(gt, point, 2023, p$size)),
                     terra::values(ref))
  }
  # A point in the zone's CRS lands on the same grid without a round trip.
  utm <- terra::project(forms[[3]], "EPSG:32631")
  r <- tessera_patch(gt, utm, 2023, p$size)
  expect_equal(as.vector(terra::ext(r)), as.vector(terra::ext(ref)))

  expect_error(tessera_patch(gt, rbind(c(0, 52), c(0.01, 52)), 2023, 8), "single point")
  expect_error(tessera_patch(gt, c(NA, 52), 2023, 8), "WGS84 bounds")
})

test_that("a seam patch follows the zone of the point's CRS", {
  gt <- fixture_store()
  p <- golden()$patch_seam
  v <- terra::vect(cbind(p$lon, p$lat), crs = "EPSG:4326")
  r31 <- tessera_patch(gt, v, 2023, p$size)
  r30 <- tessera_patch(gt, terra::project(v, "EPSG:32630"), 2023, p$size)
  expect_true(same_crs(terra::crs(r31), "EPSG:32631"))
  expect_true(same_crs(terra::crs(r30), "EPSG:32630"))
  expect_false(is.null(grid_offsets(zone_open(gt, 30), terra::ext(r30), "EPSG:32630")))
  # dst_crs takes precedence over the point's CRS.
  r <- tessera_patch(gt, terra::project(v, "EPSG:32630"), 2023, p$size,
                     dst_crs = "EPSG:32631")
  expect_identical(terra::values(r), terra::values(r31))
})

test_that("batched samples match the Python reader", {
  withr::local_options(geotessera.batch_points = 3L)
  gt <- fixture_store()
  g <- golden()
  pts <- golden_points(g)
  expect_gt(nrow(pts), 3L)
  batches <- point_batches(pts, 3L)
  expect_true(all(lengths(batches) <= 3L))
  expect_identical(sort(unlist(batches)), seq_len(nrow(pts)))
  res <- tessera_probe(gt, pts, 2024)
  expect_identical(as.character(res$status), g$probe_2024)
  expect_values(band_matrix(res), g$sample_2024)
  withr::local_options(geotessera.batch_points = 0)
  expect_error(tessera_sample(gt, pts, 2024), "geotessera.batch_points")
})

test_that("batches run on mirai daemons", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  skip_if_not(testthat::is_checking(), "the installed package may be stale")
  withr::local_options(geotessera.batch_points = 3L)
  mirai::daemons(2)
  withr::defer(mirai::daemons(0))
  gt <- fixture_store()
  g <- golden()
  res <- tessera_probe(gt, golden_points(g), 2024)
  expect_identical(as.character(res$status), g$probe_2024)
  expect_values(band_matrix(res), g$sample_2024)

  dir <- withr::local_tempdir()
  file.copy(tessera_example(), dir, recursive = TRUE)
  gone <- tessera_open(file.path(dir, "example.zarr"))
  unlink(file.path(dir, "example.zarr"), recursive = TRUE)
  expect_error(tessera_sample(gone, golden_points(g), 2024),
               class = "geotessera_error")

  fresh <- tessera_open(tessera_example())
  withr::local_options(geotessera.cache_bytes = "bad")
  expect_error(tessera_sample(fresh, golden_points(g), 2024), "max_size must be a number")
})
