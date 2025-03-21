# Now

- compressor
    - catch malloc/realloc errors
    - return object with compress fn, output stats, rankings, etc
    - allow persist/load
    - if density > 50%, flip bits before/after processing (so indexing works on
      smallest number)
    - multithreaded/configurable encode and decode
    - sanitizer/valgrind
    - can we replace malloc with lua_newuserdata across the board?
    - any opportunities to eliminate extra allocations?
    - are p_y_given_x/etc actually storing (x) and (1 - x) separately when this
      is not strictly necessary?
    - migrate some of these enhancements to tsetlin

- bitmap.create from raw string
- matrix/raw_matrix accept bitmaps or strings

# Eventually

- Consider converting this to be entirely "matrix" based

- Allow specifying int type to use for raw bitmaps
