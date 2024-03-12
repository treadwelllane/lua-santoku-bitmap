local bm = require("santoku.bitmap.capi")
local create = bm.create
local extend = bm.extend
local chunk_bits = bm.chunk_bits

local tbl = require("santoku.table")
local merge = tbl.merge

local function raw_matrix (bs, n, s, e)
  local b0 = create()
  local step = n + chunk_bits - (n % chunk_bits)
  local idx = 1
  s = s or 1
  e = e or #bs
  for i = s, e do
    extend(b0, bs[i], idx)
    idx = idx + step
  end
  return bm.raw(b0, #bs * step)
end

return merge({
  raw_matrix = raw_matrix
}, bm)
