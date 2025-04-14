# Now

- Add early stopping logic to C, normalized by mean:
    `fabs(avg(tcs, -10, -4) - avg(tcs, -5, -1)) / avg(tcs, -5, -1) < eps_norm`

- `bits` api that deals entirely with "sorted list of set bits" representation,
  using kvec, with conversion to/from roaring, packed/raw, and bits
  representations

- Bitmap: threaded top_mi calculation
- Compressor: threaded tiling of set bits
- Compressor: handle the case where n_hidden is too small, resulting in NaNs
- Compressor: compress should return list of set bits representation, threaded
- Integrate better error messages from tsetlin

# Later

- Consider switching to floats, but avoid NaNs by clamping log probabilities
  within +/- 50  or 80 (make this a compilation flag)
- Check errors for various pthread calls
- Ensure that shrink does as much as possible
- Can we eliminate one of the *2's for pc and lm? What about log_pyx_unnorm?
- Can we eliminate the pc buffer entirely?

- if density > 50%, flip bits before/after processing (so indexing works on
  smallest number)
- bitmap.create from raw string
- matrix/raw_matrix accept bitmaps or strings

# Eventually

- Clean up handling of signed/unsigned/etc integers (there are likely bugs when > UINT32_MAX)

- Consider migrating compressor to tsetlin library
- Consider converting this to be entirely "matrix" based
- Allow specifying int type to use for raw bitmaps
