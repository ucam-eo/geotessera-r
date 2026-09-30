# Byte-range fetches over HTTP and from local stores, and the chunk cache.

range_header <- function(r) {
  if (length(r) == 1L) {
    if (r >= 0) sprintf("bytes=%.0f-", r) else sprintf("bytes=-%.0f", -r)
  } else {
    sprintf("bytes=%.0f-%.0f", r[1L], r[2L] - 1)
  }
}

# A server may ignore the Range header and send the whole object.
slice_range <- function(body, r) {
  if (is.null(r)) return(body)
  n <- length(body)
  if (length(r) == 1L) {
    from <- if (r >= 0) r else max(0, n + r)
    to <- n
  } else {
    from <- r[1L]
    to <- min(r[2L], n)
  }
  if (to <= from) raw(0) else body[(from + 1):to]
}

# Fetch `urls[i]` restricted to `ranges[[i]]` concurrently. `NULL` ranges
# fetch whole objects. Missing objects return `NULL`. Stalled or failed
# requests are retried with backoff; httr2 retries only HTTP statuses of
# parallel requests, not transport failures.
http_get_ranges <- function(urls, ranges) {
  retries <- getOption("geotessera.retries", 5L)
  agent <- paste0("geotessera-r/", utils::packageVersion("geotessera"))
  out <- vector("list", length(urls))
  pending <- seq_along(urls)
  for (attempt in 0:retries) {
    if (attempt) Sys.sleep(min(30, 2^(attempt - 1)) * retry_jitter())
    reqs <- lapply(pending, function(i) {
      req <- httr2::request(urls[i])
      req <- httr2::req_user_agent(req, agent)
      req <- httr2::req_options(req, timeout = 120, connecttimeout = 30,
                                low_speed_limit = 16384, low_speed_time = 10)
      req <- httr2::req_error(req, is_error = function(resp) FALSE)
      r <- ranges[[i]]
      if (is.null(r)) req else httr2::req_headers(req, Range = range_header(r))
    })
    resps <- httr2::req_perform_parallel(
      reqs, on_error = "continue", progress = FALSE,
      max_active = getOption("geotessera.concurrency", 32L)
    )
    failed <- logical(length(pending))
    why <- character(length(pending))
    for (j in seq_along(pending)) {
      i <- pending[j]
      resp <- resps[[j]]
      if (inherits(resp, "error")) {
        failed[j] <- TRUE
        why[j] <- conditionMessage(resp)
        next
      }
      status <- httr2::resp_status(resp)
      if (status == 200L) {
        out[i] <- list(slice_range(httr2::resp_body_raw(resp), ranges[[i]]))
      } else if (status == 206L) {
        out[i] <- list(httr2::resp_body_raw(resp))
      } else if (status %in% c(408L, 425L, 429L) || status >= 500L) {
        failed[j] <- TRUE
        why[j] <- paste("HTTP", status)
      } else if (status != 404L) {
        abort(sprintf("HTTP %d fetching %s", status, urls[i]),
              class = "geotessera_http_error")
      }
    }
    if (!any(failed)) return(out)
    first <- which(failed)[1L]
    pending <- pending[failed]
  }
  abort(sprintf("Failed to fetch %s after %d retries: %s",
                urls[pending[1L]], retries, why[first]),
        class = "geotessera_http_error")
}

# Fetch byte ranges of store keys, concurrently where the store allows.
store_get_ranges <- function(store, keys, ranges) {
  if (inherits(store, "tessera_http_store")) {
    store$get_ranges(keys, ranges)
  } else {
    Map(function(key, r) store$get(key, byte_range = r), keys, ranges)
  }
}

# A size-bounded LRU cache of shard indices and compressed chunks, in
# memory and, given `dir`, backed by disk so that downloads persist between
# sessions. cachem restricts key characters, so keys are hashed.
new_cache <- function(max_bytes, dir = NULL) {
  mem <- cachem::cache_mem(max_size = max_bytes, evict = "lru")
  if (is.null(dir)) return(mem)
  disk <- cachem::cache_disk(
    dir, evict = "lru",
    max_size = getOption("geotessera.disk_cache_bytes", 10 * 1024^3)
  )
  cachem::cache_layered(mem, disk)
}

cache_get <- function(cache, key) {
  value <- cache$get(digest::digest(key, "xxhash64", serialize = FALSE))
  if (cachem::is.key_missing(value)) NULL else value
}

cache_put <- function(cache, key, value) {
  cache$set(digest::digest(key, "xxhash64", serialize = FALSE), value)
}

# A factor in [0.5, 1) that spreads out retries without advancing the
# user's random number stream.
retry_jitter <- function() 0.5 + (as.numeric(Sys.time()) * 1000) %% 1000 / 2000
