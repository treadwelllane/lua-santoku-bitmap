# Now

- Separate create and train functions
- Rename compress to encode
- Determine stopping condition in callback, not in C

# Later

- are p_y_given_x/etc actually storing (x) and (1 - x) separately when this
  is not strictly necessary?
- if density > 50%, flip bits before/after processing (so indexing works on
  smallest number)
- fn to get stats, feature rankings, etc
- bitmap.create from raw string
- matrix/raw_matrix accept bitmaps or strings

# Eventually

- Clean up handling of signed/unsigned/etc integers (there are likely bugs when > UINT32_MAX)

- Consider migrating compressor to tsetlin library
- Consider converting this to be entirely "matrix" based
- Allow specifying int type to use for raw bitmaps
