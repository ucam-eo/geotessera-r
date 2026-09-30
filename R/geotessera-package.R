#' geotessera: read Tessera embeddings from Zarr stores
#'
#' Tessera is a foundation model for Earth observation that summarises a
#' year of Sentinel-1 and Sentinel-2 observations as an embedding at each
#' 10 m pixel. This package reads published embeddings from Zarr stores.
#'
#' Open a store with [tessera_open()], by default the v1.1 global dClimate
#' run; [tessera_datasets()] lists the others. Then:
#' * sample points with [tessera_sample()] or [tessera_probe()];
#' * read a region with [tessera_region()] or [tessera_map_region()];
#' * read a fixed-size patch with [tessera_patch()];
#' * write GeoTIFF files with [tessera_export_geotiffs()].
#'
#' Coordinates are WGS84 longitude and latitude unless a spatial object
#' declares another CRS. Embeddings are returned on the native UTM grid of
#' the zone that holds them and are never resampled, except for the parts of
#' a [tessera_patch()] that cross into another zone or onto a `dst_crs`.
#' Band indices are one-based.
#'
#' Errors signalled by the package have class `geotessera_error`, and failed
#' downloads also `geotessera_http_error`.
#'
#' @section Options:
#' \describe{
#'   \item{`geotessera.cache_dir`}{The default `cache_dir` of
#'     [tessera_open()], where remote stores keep metadata and downloaded
#'     chunks between sessions. Unset by default.}
#'   \item{`geotessera.disk_cache_bytes`}{Size bound in bytes of each
#'     store's downloaded chunks in `cache_dir`, least recently used first
#'     out. The default is 10 GiB.}
#'   \item{`geotessera.concurrency`}{Maximum concurrent HTTP requests.
#'     The default is 32.}
#'   \item{`geotessera.retries`}{Retries of a failed or stalled HTTP
#'     request, with exponential backoff. The default is 5.}
#'   \item{`geotessera.cache_bytes`}{Size bound in bytes of the in-memory
#'     cache of shard indices and compressed chunks, set when a store is
#'     opened. The default is 256 MiB.}
#' }
#'
#' @references Feng Z. et al. (2025). TESSERA: Temporal Embeddings of
#'   Surface Spectra for Earth Representation and Analysis.
#'   \doi{10.48550/arXiv.2506.20380}
#'
#' @keywords internal
#' @importFrom R6 R6Class
#' @importFrom curl new_handle
#' @importFrom digest digest
"_PACKAGE"

# curl and digest are imported so that zarr, which loads them with a
# message, finds them already loaded.
