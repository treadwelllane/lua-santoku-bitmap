local bm = require("santoku.bitmap.capi")
local create = bm.create
local extend = bm.extend

local tbl = require("santoku.table")
local merge = tbl.merge

local function raw_matrix (bs, step, s, e)
  local b0 = create()
  local idx = 1
  s = s or 1
  e = e or #bs
  for i = s, e do
    extend(b0, bs[i], idx)
    idx = idx + step
  end
  return bm.raw(b0, #bs * step)
end

local function next (b, i)
  for j = i + 1, bm.maximum(b) do
    if bm.get(b, j) then
      return j
    end
  end
end

local function hamming (a, b, t)
  if t then
    bm.clear(t)
  else
    t = bm.create()
  end
  bm["or"](t, a)
  bm["xor"](t, b)
  return bm.cardinality(t)
end

return merge({
  raw_matrix = raw_matrix,
  next = next,
  hamming = hamming,
}, bm)
