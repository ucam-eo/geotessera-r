fixture_store <- function() {
  tessera_open(tessera_example())
}

# The band columns added by tessera_sample() and tessera_probe().
band_matrix <- function(df) {
  as.matrix(as.data.frame(df)[grep("^band_", names(df))])
}

golden <- function() {
  jsonlite::fromJSON(test_path("fixtures", "golden.json"))
}

golden_points <- function(g) {
  unname(do.call(rbind, g$points))
}

# Python returns (row, column, band) float32 arrays; values agree to float32
# precision.
expect_values <- function(object, expected, tolerance = 1e-5) {
  object <- unname(object)
  expected <- unname(expected)
  expect_equal(dim(object), dim(expected))
  expect_identical(is.na(object), is.na(expected))
  expect_equal(object[!is.na(object)], expected[!is.na(expected)],
               tolerance = tolerance)
}
