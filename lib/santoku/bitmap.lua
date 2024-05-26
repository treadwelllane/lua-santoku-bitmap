local bm = require("santoku.bitmap.capi")
local num = require("santoku.num")
local tbl = require("santoku.table")

local function matrix (bs, step, s, e, rows)
  local b0 = bm.create()
  step = num.ceil(step / bm.chunk_bits) * bm.chunk_bits
  local idx = 1
  s = s or 1
  e = e or #bs
  local n = 0
  for i = s, e do
    bm.extend(b0, bs[i], idx)
    idx = idx + step
    n = n + 1
    if rows and n >= rows then
      break
    end
  end
  return b0, step
end

local function raw_matrix (bs, step0, s, e, rows)
  local b0, step1 = matrix(bs, step0, s, e, rows)
  local l = rows or #bs
  return bm.raw(b0, l * step1)
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

return tbl.merge({
  raw_matrix = raw_matrix,
  matrix = matrix,
  next = next,
  hamming = hamming,
}, bm)
