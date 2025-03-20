local bm = require("santoku.bitmap.capi")
local arr = require("santoku.array")
local tbl = require("santoku.table")

local function matrix (bs, step, s, e, rows)
  local b0 = bm.create()
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

local function tostring (b, n)
  n = n or bm.cardinality(b)
  n = n % 32 == 0 and n or (n + (32 - (n % 32)))
  local x = 1
  local t = bm.bits(b, n)
  local out = {}
  for i = 1, #t do
    if i > n then
      break
    end
    local bit = t[i]
    for j = x, bit - 1 do
      out[j] = "0"
    end
    out[bit] = "1"
    x = bit + 1
  end
  while x <= n do
    out[x] = "0"
    x = x + 1
  end
  return arr.concat(out)
end

return tbl.merge({
  tostring = tostring,
  raw_matrix = raw_matrix,
  matrix = matrix,
  next = next,
  hamming = hamming,
}, bm)
