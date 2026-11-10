# geotessera

geotessera reads [Tessera](https://github.com/ucam-eo/tessera) embeddings
from local or remote Zarr stores in R. Tessera summarises a year of
Sentinel-1 and Sentinel-2 observations as an embedding at each 10 m pixel.

Queries take WGS84 longitudes and latitudes, or `sf` and `terra` objects in
any CRS. Embeddings are dequantized and returned as columns added to your
points, or as `terra` rasters on the native UTM grid of the zone that holds
them. Band indices are one-based.

See `vignette("geotessera")` for a tour that runs offline on the example
store from `tessera_example()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ucam-eo/geotessera-r")
```

The package requires R 4.4 or later and uses
[zarr](https://cran.r-project.org/package=zarr) to decode arrays and
[terra](https://cran.r-project.org/package=terra) for rasters and
projections.

## Open a store

```r
library(geotessera)

gt <- tessera_open()
gt$years
```

`tessera_open()` reads the public v1.1 store of the global dClimate run.
`tessera_datasets()` lists the published datasets; select another with
`tessera_store_url()`, or pass the path of a local store:

```r
tessera_datasets()
gt1 <- tessera_open(tessera_store_url("v1"))
cam <- tessera_open(tessera_store_url("v1.1", "cambridge"))
gt2 <- tessera_open(tessera_store_url("v2"))
local <- tessera_open("/data/tessera.zarr")
```

Each dataset is a separate inference run, a version and a variant.
Embeddings of different datasets are unrelated, even within one version,
so use one dataset per analysis: the v1.1 `dclimate` and `cambridge`
variants share a model but not a run. Stores must carry consolidated Zarr
metadata.

Set `cache_dir`, or the `geotessera.cache_dir` option, to keep store
metadata and downloaded chunks on disk between sessions, so repeated reads
of an area are not downloaded again:

```r
options(geotessera.cache_dir = tools::R_user_dir("geotessera", "cache"))
```

Each store URL has its own subdirectory, with chunks bounded by the
`geotessera.disk_cache_bytes` option (default 10 GiB). Cached data is never
revalidated; delete a store's subdirectory to pick up changes to the
store.

## Sample points

```r
pts <- data.frame(site = c("Cambridge", "Liverpool"),
                  lon = c(0.12, -2.97), lat = c(52.20, 53.44))
tessera_sample(gt, pts, 2024)
```

`tessera_sample()` returns the points with columns `band_1`, `band_2` and
so on added, `NA` for points without an embedding. Data frames, tibbles,
`sf` and `SpatVector` point layers keep their class; a matrix of longitudes
and latitudes becomes a data frame.

`tessera_probe()` also adds a `status` factor: `valid`, `water`, `nodata`
for an unwritten pixel, or `outside` beyond coverage.

```r
table(tessera_probe(gt, pts, 2024)$status)
```

### Many points

Embeddings are stored in chunks of 32 by 32 pixels and every band, about
110 KB each, and a point is read by downloading its whole chunk. Points
within a few hundred metres of each other share chunks; a million scattered
points download on the order of 100 GB, so large jobs are bound by
bandwidth rather than CPU.

Points are read in batches of nearby points, which bounds memory use. With
[mirai](https://mirai.r-lib.org) daemons set, batches are read in parallel:

```r
mirai::daemons(4)
gt <- tessera_open(cache_dir = tools::R_user_dir("geotessera", "cache"))
res <- tessera_sample(gt, pts, 2024)
mirai::daemons(0)
```

Each daemon opens the store itself and needs geotessera installed. A few
daemons usually saturate the connection. Give `tessera_open()` a
`cache_dir` so that daemons and later runs share downloads, and raise
`geotessera.disk_cache_bytes` above the size of the job. Where points are
dense, reading the area with `tessera_region()` and extracting with
`terra::extract()` can download less. Running close to the store, on a
cloud machine near its storage, helps more than more daemons.

## Read a region

```r
bbox <- c(0.05, 52.15, 0.20, 52.25)
r <- tessera_region(gt, bbox, 2024)
terra::plotRGB(terra::stretch(r[[1:3]]))
```

`tessera_region()` takes bounds as west, south, east and north, or as an
`sf` or `terra` object whose extent is used. It returns a `SpatRaster` on
the grid of one UTM zone: the zone whose CRS the bounds are given in, or
else the zone holding their centre.
`tessera_region_quantized()` returns the stored integers with a `scale`
layer, which `tessera_dequantize()` converts.

`tessera_map_region()` applies a function to the region in row strips of
at most `strip_rows` rows, so the region need not fit in memory:

```r
means <- tessera_map_region(gt, bbox, 2024, terra::global, "mean", na.rm = TRUE)
```

## Read a patch

```r
p <- tessera_patch(gt, c(0.12, 52.20), 2024, size_px = 256)
```

The point, given in any form accepted by `tessera_sample()`, falls in pixel
`[size_px %/% 2 + 1, size_px %/% 2 + 1]`. The patch lies on the native grid
of the point's UTM zone, or of the zone whose CRS the point is given in. Where it crosses into
a neighbouring zone, only the pixels beyond the seam are resampled, by
nearest neighbour. Set `dst_crs` to lay the patch on another projected CRS.

## Working on the native grid

Tessera embeddings are stored on the 10 m UTM grid of each zone, and the
package returns them on that grid without resampling. Bounds and points
given in a zone's CRS are used as given; those in other CRSs are projected
once, straight to the zone's CRS. To read exactly the grid of an earlier
result, pass it as the bounds:

```r
r <- tessera_region(gt, bbox, 2024)
r2023 <- tessera_region(gt, r, 2023)   # identical grid to r
```

Reproject your own data to the embeddings' grid, rather than the reverse,
for example with `terra::project(other, r, method = "near")`.

## Export GeoTIFFs

```r
tessera_export_geotiffs(gt, bbox, 2024, "region", dry_run = TRUE)
files <- tessera_export_geotiffs(gt, bbox, 2024, "region", bands = 1:3)
```

Each intersecting UTM zone is written to `tessera_YEAR_utmNN.tif` as
float32 with NaN nodata and the same `band_N` layer names as
`tessera_region()`.

Exports of a published dataset record it in `TESSERA_DATASET_*` tags on
each file and in a `tessera_metadata.json` file in the output directory,
as the Python package does. Exporting another dataset into that directory
is an error, so embeddings of different runs are never mixed.

## Matryoshka depths

v2 stores publish prefixes of the full embedding. Select one with `depth`:

```r
x16 <- tessera_sample(gt2, pts, 2024, depth = 16)
```

## Options

`geotessera.concurrency` sets the number of concurrent HTTP requests
(default 32). `geotessera.retries` sets the retries for a failed request
(default 5). `geotessera.cache_dir` sets the default `cache_dir` of
`tessera_open()`, and `geotessera.disk_cache_bytes` bounds each store's
chunks there (default 10 GiB). `geotessera.cache_bytes` bounds the in-memory cache of shard
indices and compressed chunks (default 256 MiB).
`geotessera.batch_points` sets the largest batch of points read at once
(default 10000).

## Tests

```sh
Rscript -e 'testthat::test_local()'
GEOTESSERA_ONLINE=true Rscript -e 'testthat::test_local()'
```

The offline tests compare results with the Python reader on the synthetic
example store in `inst/extdata`. Regenerate the store and expected values
from a `geotessera` checkout:

```sh
uv run --project ../geotessera python data-raw/make-fixture.py
```

## Citation

If you use Tessera embeddings in research, cite the model paper;
`citation("geotessera")` gives the reference.

## Licence

MIT. See `LICENSE.md`.
