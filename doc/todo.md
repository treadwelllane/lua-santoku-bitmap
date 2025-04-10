# Now

- Handle the case where n_hidden is too small, resulting in NaNs
- Potential memory leak when pinning nodes

# Later

- Consider switching to floats, but avoid NaNs by clamping log probabilities
  within +/- 50  or 80 (make this a compilation flag)
- Check errors for various pthread calls
- Ensure that shrink does as much as possible
- Can we eliminate one of the *2's for pc and lm? What about log_pyx_unnorm?
- Can we eliminate the pc buffer entirely?

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
