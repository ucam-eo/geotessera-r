# geotessera 0.2.0

* Requires zarr 1.0.0+ and R 4.4+. `tessera_open()` failed with zarr 1.0 (#1).
* Sampling many points is faster as planning byte ranges no longer takes
  time quadratic in the number of chunks (reported by Jayden Engbert). 
* `tessera_sample()` and `tessera_probe()` read points in batches of nearby
  points (`geotessera.batch_points`, default 10000).
  Batches run in parallel when mirai daemons are set with `mirai::daemons()`.

# geotessera 0.1.0

* Initial release. Reads Tessera embeddings from local and remote Zarr v3
  stores, a port of the `GeoTesseraZarr` reader of the Python `geotessera`
  package.
* The default store is the v1.1 global dClimate run. `tessera_datasets()`
  lists the published datasets and `tessera_store_url()` selects one by
  version and variant.
* Exports of a published dataset record it in GeoTIFF tags and a
  `tessera_metadata.json` file, and refuse directories holding another
  dataset.
* `cache_dir` keeps downloaded chunks on disk as well as store metadata,
  bounded by the `geotessera.disk_cache_bytes` option.
* Reads warn when a hemisphere of a zone has no complete embeddings for
  the requested year.
