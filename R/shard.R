# Reads of sharded Zarr arrays. Array metadata and the codec pipelines come
# from the zarr package; this layer caches shard indices and compressed
# chunks, and fetches the byte ranges of many inner chunks concurrently.
#
# zarr_array$read() is not used because, as of zarr 0.5.1, its sharded
# reader keeps no state between calls: every read downloads the shard index
# again (256 KiB for a 4096-pixel shard of 32-pixel chunks) and decodes its
# chunks again, as Zarr.options$chunk_cache_bytes applies only to unsharded
# arrays. It also fetches shards one at a time, and one byte span from the
# first to the last inner chunk needed, which for a narrow window crosses
# most of the shard's width. Against source.coop, sampling 200 points took
# 536 s rather than 19 s, and reading a 1024 by 512 pixel window in four
# row strips downloaded 675 MB rather than 61 MB.
#
# Indices are 0-based and follow the Zarr dimension order. Decoded chunks and
# assembled blocks are R arrays with the dimensions reversed, so an element
# of a (time, band, y, x) array is addressed as [x, y, band, time].

sharded_reader <- function(arr, store, cache) {
  meta <- arr$metadata
  sharding <- meta$codecs[[1L]]
  if (length(meta$codecs) != 1L || !identical(sharding$name, "sharding_indexed")) {
    abort(sprintf("Array '%s' is not sharded", arr$path))
  }
  chunking <- arr$chunking
  codecs <- chunking$codecs
  if (!identical(codecs[[1L]]$name, "bytes")) {
    abort(sprintf("Unsupported codec '%s' in array '%s'", codecs[[1L]]$name, arr$path))
  }
  inner <- as.numeric(chunking$inner_shape)
  shard <- as.numeric(chunking$chunk_shape)
  grid <- shard / inner
  index_codecs <- chunking$index_codecs
  crc <- any(vapply(index_codecs, function(c) c$name == "crc32c", logical(1)))
  encoding <- chunking$chunk_encoding
  list(
    path = arr$path,
    shape = as.numeric(meta$shape),
    inner = inner,
    grid = grid,
    strides = c(rev(cumprod(rev(grid)))[-1L], 1),
    n_inner = prod(grid),
    index_bytes = prod(grid) * 16 + if (crc) 4 else 0,
    index_at_end = !identical(sharding$configuration$index_location, "start"),
    codecs = codecs,
    index_codecs = index_codecs,
    key_prefix = paste0(arr$prefix, encoding$pre),
    sep = encoding$sep,
    fill = arr$data_type$fill_value,
    store = store,
    cache = cache
  )
}

decode_shard_index <- function(reader, bytes) {
  for (codec in rev(reader$index_codecs)) {
    if (codec$name != "bytes") bytes <- codec$decode(bytes)
  }
  if (length(bytes) != reader$n_inner * 16) {
    abort(sprintf("Malformed shard index in array '%s'", reader$path))
  }
  w <- matrix(readBin(bytes, "integer", n = length(bytes) / 2, size = 2,
                      signed = FALSE, endian = "little"), nrow = 8L)
  scale <- c(1, 65536, 65536^2, 65536^3)
  offset <- colSums(w[1:4, , drop = FALSE] * scale)
  length <- colSums(w[5:8, , drop = FALSE] * scale)
  offset[w[4L, ] == 65535 & w[8L, ] == 65535] <- NA
  list(offset = offset, length = length)
}

# Store keys of the shards holding inner chunks, and each chunk's 1-based
# position in its shard index.
shard_location <- function(reader, cc) {
  cc <- matrix(cc, ncol = length(reader$inner))
  grid <- matrix(reader$grid, nrow(cc), ncol(cc), byrow = TRUE)
  shard_cc <- floor(cc / grid)
  linear <- as.integer(as.vector((cc - shard_cc * grid) %*% reader$strides) + 1)
  storage.mode(shard_cc) <- "integer"
  keys <- paste0(reader$key_prefix,
                 do.call(paste, c(as.data.frame(shard_cc), sep = reader$sep)))
  list(keys = keys, linear = linear)
}

# Load uncached shard indices of several readers in one concurrent round.
# `wanted` pairs each reader with the store keys of its shards. A missing
# shard has index `FALSE`.
load_indices <- function(wanted) {
  todo <- list()
  for (w in wanted) {
    for (k in unique(w$keys)) {
      if (is.null(cache_get(w$reader$cache, paste0("index:", k)))) {
        todo[[length(todo) + 1L]] <- list(reader = w$reader, key = k)
      }
    }
  }
  if (!length(todo)) return(invisible())
  ranges <- lapply(todo, function(t) {
    if (t$reader$index_at_end) -t$reader$index_bytes else c(0, t$reader$index_bytes)
  })
  keys <- vapply(todo, `[[`, character(1), "key")
  bodies <- store_get_ranges(todo[[1L]]$reader$store, keys, ranges)
  for (j in seq_along(todo)) {
    reader <- todo[[j]]$reader
    index <- if (is.null(bodies[[j]])) FALSE else decode_shard_index(reader, bodies[[j]])
    cache_put(reader$cache, paste0("index:", keys[j]), index)
  }
  invisible()
}

# Compressed bytes of inner chunks, one row of `cc` per chunk in inner-chunk
# grid coordinates. Absent chunks are `NULL`.
fetch_chunks <- function(reader, cc) {
  cc <- matrix(cc, ncol = length(reader$inner))
  n <- nrow(cc)
  out <- vector("list", n)
  if (!n) return(out)
  loc <- shard_location(reader, cc)
  keys <- loc$keys
  linear <- loc$linear
  load_indices(list(list(reader = reader, keys = keys)))

  chunk_keys <- paste0(keys, "#", linear)
  offset <- rep(NA_real_, n)
  length <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    cached <- cache_get(reader$cache, chunk_keys[i])
    if (!is.null(cached)) {
      out[i] <- list(cached)
      next
    }
    index <- cache_get(reader$cache, paste0("index:", keys[i]))
    if (isFALSE(index) || is.na(index$offset[linear[i]])) next
    offset[i] <- index$offset[linear[i]]
    length[i] <- index$length[linear[i]]
  }

  todo <- which(!is.na(offset))
  if (!length(todo)) return(out)
  spans <- plan_spans(keys[todo], offset[todo], length[todo])
  bodies <- store_get_ranges(reader$store, spans$key,
                             Map(c, spans$start, spans$end))
  for (s in seq_len(nrow(spans))) {
    body <- bodies[[s]]
    if (is.null(body)) {
      abort(sprintf("Shard '%s' disappeared during a read", spans$key[s]))
    }
    for (i in todo[spans$members[[s]]]) {
      from <- offset[i] - spans$start[s]
      bytes <- body[(from + 1):(from + length[i])]
      cache_put(reader$cache, chunk_keys[i], bytes)
      out[i] <- list(bytes)
    }
  }
  out
}

# Coalesce nearby byte ranges of the same shard into single requests. Public
# object stores can stall on large ranges, so spans stay small.
plan_spans <- function(keys, offset, length, max_gap = 16 * 1024,
                       max_span = 1024^2) {
  ord <- order(keys, offset)
  key <- character(0)
  start <- numeric(0)
  end <- numeric(0)
  members <- list()
  for (i in ord) {
    s <- length(key)
    if (s && key[s] == keys[i] && offset[i] - end[s] <= max_gap &&
        offset[i] + length[i] - start[s] <= max_span) {
      end[s] <- max(end[s], offset[i] + length[i])
      members[[s]] <- c(members[[s]], i)
    } else {
      key <- c(key, keys[i])
      start <- c(start, offset[i])
      end <- c(end, offset[i] + length[i])
      members[[s + 1L]] <- i
    }
  }
  spans <- data.frame(key = key, start = start, end = end)
  spans$members <- members
  spans
}

decode_chunk <- function(reader, bytes) {
  if (is.null(bytes)) {
    return(array(reader$fill, dim = rev(reader$inner)))
  }
  for (codec in rev(reader$codecs)) bytes <- codec$decode(bytes)
  dim(bytes) <- rev(reader$inner)
  bytes
}

# The block [lo, hi) as an array with reversed dimensions.
read_block <- function(reader, lo, hi) {
  nd <- length(lo)
  c0 <- floor(lo / reader$inner)
  c1 <- floor((hi - 1) / reader$inner)
  cc <- as.matrix(expand.grid(lapply(seq_len(nd), function(d) c0[d]:c1[d])))
  chunks <- fetch_chunks(reader, cc)
  out <- array(reader$fill, dim = rev(hi - lo))
  for (i in seq_len(nrow(cc))) {
    if (is.null(chunks[[i]])) next
    buf <- decode_chunk(reader, chunks[[i]])
    origin <- cc[i, ] * reader$inner
    s <- pmax(lo, origin)
    e <- pmin(hi, origin + reader$inner)
    src <- lapply(nd:1, function(d) (s[d] - origin[d] + 1):(e[d] - origin[d]))
    dst <- lapply(nd:1, function(d) (s[d] - lo[d] + 1):(e[d] - lo[d]))
    # Assigned in place here; passing `out` to a helper would copy it.
    switch(
      nd,
      out[dst[[1]]] <- buf[src[[1]]],
      out[dst[[1]], dst[[2]]] <- buf[src[[1]], src[[2]]],
      out[dst[[1]], dst[[2]], dst[[3]]] <- buf[src[[1]], src[[2]], src[[3]]],
      out[dst[[1]], dst[[2]], dst[[3]], dst[[4]]] <-
        buf[src[[1]], src[[2]], src[[3]], src[[4]]]
    )
  }
  out
}

# Values at points, one row of `idx` (0-based, Zarr order) per point. Dims in
# `whole` are read in full and must lie within one inner chunk; the result
# has one row per point and one column per element along them.
read_points <- function(reader, idx, whole = integer(0)) {
  nd <- length(reader$inner)
  idx <- matrix(idx, ncol = nd)
  n <- nrow(idx)
  width <- prod(reader$shape[whole])
  out <- matrix(if (is.integer(reader$fill)) NA_integer_ else NA_real_, n, width)
  if (!n) return(out)
  if (any(reader$inner[whole] < reader$shape[whole])) {
    abort(sprintf("Array '%s' splits a whole dimension across chunks", reader$path))
  }
  idx[, whole] <- 0
  cc <- floor(idx / matrix(reader$inner, n, nd, byrow = TRUE))
  id <- as.vector(cc %*% c(rev(cumprod(rev(ceiling(reader$shape / reader$inner))))[-1L], 1))
  groups <- split(seq_len(n), id)
  first <- vapply(groups, `[`, integer(1), 1L)
  chunks <- fetch_chunks(reader, cc[first, , drop = FALSE])
  full <- as.matrix(expand.grid(lapply(whole, function(d) seq_len(reader$shape[d]))))
  n_full <- if (length(whole)) nrow(full) else 1L
  for (g in seq_along(groups)) {
    members <- groups[[g]]
    buf <- decode_chunk(reader, chunks[[g]])
    local <- idx[members, , drop = FALSE] -
      matrix(cc[first[g], ] * reader$inner, length(members), nd, byrow = TRUE) + 1
    sel <- local[rep(seq_along(members), times = n_full), , drop = FALSE]
    if (length(whole)) sel[, whole] <- full[rep(seq_len(n_full), each = length(members)), ]
    out[members, ] <- buf[sel[, nd:1, drop = FALSE]]
  }
  out
}
