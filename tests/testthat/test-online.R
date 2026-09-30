# Reads the public store; set GEOTESSERA_ONLINE=true to run. Expected values
# come from the Python reader.

test_that("the public v1 store matches the Python reader", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("GEOTESSERA_ONLINE"), "true"),
              "GEOTESSERA_ONLINE is not true")
  skip_if_offline("data.source.coop")

  cache_dir <- withr::local_tempfile()
  gt <- tessera_open(tessera_store_url("v1"), cache_dir = cache_dir)
  expect_identical(gt$dataset, list(version = "1.0", variant = "vultr"))
  expect_identical(gt$years, 2017:2025)
  expect_identical(gt$n_bands, 128L)

  pts <- rbind(c(0.12, 52.20), c(-2.97, 53.44), c(-10.0, 55.5))
  res <- tessera_probe(gt, pts, 2024)
  expect_identical(as.character(res$status), c("valid", "valid", "nodata"))
  emb <- band_matrix(res)
  expect_equal(unname(emb[1, 1:3]),
               c(8.03232479095459, -0.3162332773208618, 0.505973219871521),
               tolerance = 1e-6)
  expect_equal(unname(emb[2, 1:3]),
               c(3.368870973587036, -0.6501330137252808, -0.17730900645256042),
               tolerance = 1e-6)

  r <- tessera_region(gt, c(0.10, 52.19, 0.12, 52.20), 2024)
  expect_equal(dim(r), c(117L, 141L, 128L))
  expect_equal(as.vector(terra::ext(r))[c("xmin", "ymax")],
               c(xmin = 301780, ymax = 5787250))
  a <- terra::as.array(r)
  expect_equal(a[1, 1, 1:3], c(0.858729898929596, -1.3271280527114868, 0.31226542592048645),
               tolerance = 1e-6)
  expect_equal(mean(a), 0.18927724659442902, tolerance = 1e-6)

  # A second session reads the chunks from disk, without downloading.
  chunks <- file.path(cache_dir, "v1", "chunks")
  expect_gt(length(list.files(chunks)), 0)
  offline <- tessera_open(tessera_store_url("v1"), cache_dir = cache_dir)
  withr::local_options(geotessera.retries = 0L)
  local_mocked_bindings(http_get_ranges = function(...) stop("downloaded"))
  expect_equal(terra::as.array(tessera_region(offline, c(0.10, 52.19, 0.12, 52.20), 2024)), a)
})

test_that("the public v1.1 dClimate store matches the Python reader", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("GEOTESSERA_ONLINE"), "true"),
              "GEOTESSERA_ONLINE is not true")
  skip_if_offline("data.source.coop")

  gt <- tessera_open()
  expect_identical(gt$dataset, list(version = "1.1", variant = "dclimate"))
  expect_identical(gt$years, 2017:2025)
  expect_identical(length(gt$zone_ids), 60L)
  expect_identical(gt$model_version, "https://geotessera.org/model/1.1")

  pts <- rbind(c(0.1218, 52.2053), c(0.0953, 52.1795), c(1.9, 52.3), c(-10.0, 55.5))
  res <- tessera_probe(gt, pts, 2024)
  expect_identical(as.character(res$status), c("valid", "valid", "valid", "nodata"))
  emb <- band_matrix(res)
  expect_equal(unname(emb[1, 1:5]),
               c(1.1703, -0.2219, -0.1412, 0.2825, 0.3027), tolerance = 1e-4)
  expect_equal(unname(rowMeans(emb[1:3, ])),
               c(0.03136918321251869, -0.04314561188220978, 0.026263386011123657),
               tolerance = 1e-6)

  r <- tessera_region(gt, c(0.10, 52.19, 0.13, 52.21), 2024)
  expect_equal(dim(r), c(231L, 214L, 128L))
  expect_identical(terra::crs(r, describe = TRUE)$code, "32631")
  xy <- terra::project(pts[1, , drop = FALSE], "EPSG:4326", terra::crs(r))
  expect_equal(unname(unlist(terra::extract(r, xy)[1, ])),
               unname(emb[1, ]), tolerance = 1e-6)

  p <- tessera_patch(gt, pts[1, ], 2024, size_px = 64)
  expect_equal(terra::as.array(p)[33, 33, ], unname(emb[1, ]),
               tolerance = 1e-6)
})
