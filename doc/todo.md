# Now

- does the compress fn rely on any arrays allocated with size dependent on
  n_samples? If so, we need to realloc based on input size.

- can we replace malloc with lua_newuserdata across the board?

- compressor
    - return object with compress fn, output stats, rankings, etc
    - allow persist/load
    - index set bits per sample, which should dramatically increase perf
    - if density > 50%, flip bits before/after processing (so indexing works on
      smallest number)
    - sanitizer/valgrind
    - refactor to support auto-vectorization
    - multithreaded/configurable encode and decode
    - can we reduce number of mallocs (are some temp arrays uncessary)?

- bitmap.create from raw string
- matrix/raw_matrix accept bitmaps or strings

# Eventually

- Consider converting this to be entirely "matrix" based

- Allow specifying int type to use for raw bitmaps
