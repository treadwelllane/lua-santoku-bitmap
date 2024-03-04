local bm = require("santoku.bitmap.capi")
local create = bm.create
local extend = bm.extend
local chunk_bits = bm.chunk_bits

local it = require("santoku.iter")
local ivals = it.ivals

local tbl = require("santoku.table")
local assign = tbl.assign

local function raw_matrix (bs, n)
  local b0 = create()
  local step = n + chunk_bits - (n % chunk_bits)
  local idx = 1
  for p in ivals(bs) do
    extend(b0, p, idx)
    idx = idx + step
  end
  return bm.raw(b0, #bs * step)
end

return assign({
  raw_matrix = raw_matrix
}, bm, false)
