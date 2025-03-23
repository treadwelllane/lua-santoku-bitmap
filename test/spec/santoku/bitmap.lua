local test = require("santoku.test")
local utc = require("santoku.utc")
local err = require("santoku.error")
local assert = err.assert
local vdt = require("santoku.validate")
local eq = vdt.isequal
local bm = require("santoku.bitmap")

test("chunk bits", function ()
  assert(eq(bm.chunk_bits, 32))
end)

test("binary", function ()

  local b

  b = bm.create()
  bm.set(b, 1)
  assert(eq(bm.tostring(b), "10000000000000000000000000000000"))

  b = bm.create()
  bm.set(b, 1)
  bm.set(b, 2)
  bm.set(b, 4)
  assert(eq(bm.tostring(b), "11010000000000000000000000000000"))

  b = bm.create()
  bm.set(b, 8)
  assert(eq(bm.tostring(b), "00000001000000000000000000000000"))

  b = bm.create()
  bm.set(b, 10)
  assert(eq(bm.tostring(b), "00000000010000000000000000000000"))

end)

test("unbalanced", function ()

  local b0 = bm.create()
  bm.set(b0, 2)
  bm.set(b0, 5)
  bm.set(b0, 6)
  bm.set(b0, 21)

  assert(eq(bm.tostring(b0, 32), "01001100000000000000100000000000"))

  local b1 = bm.create()
  bm.set(b1, 2)
  bm.set(b1, 5)
  bm.set(b1, 6)
  bm.set(b1, 21)

  bm.extend(b0, b1, 33)

  assert(eq(bm.tostring(b0, 64), "0100110000000000000010000000000001001100000000000000100000000000"))

end)

test("extends", function ()

  local a = bm.create()
  bm.set(a, 1)
  bm.set(a, 4)

  local b = bm.create()
  bm.set(b, 1)
  bm.set(b, 4)

  local c = bm.create()
  bm.extend(c, a, 1)
  bm.extend(c, a, 33)

  assert(eq(bm.tostring(c, 61), "1001000000000000000000000000000010010000000000000000000000000000"))

end)

test("create from table", function ()

  local b = bm.create({ [1] = true, [3] = true, [9] = true }, 9)
  assert(eq(bm.tostring(b, 9), "10100000100000000000000000000000"))

end)

test("raw, from_raw", function ()

  local b = bm.from_raw(bm.raw(bm.create({ [1] = true, [3] = true, [9] = true }, 9)))
  assert(eq(bm.tostring(b, 9), "10100000100000000000000000000000"))

end)

test("flip", function ()

  local b = bm.create({ [1] = true, [3] = true, [9] = true }, 9)
  local c = bm.copy(b)
  bm.flip(c, 1, 9)
  assert(eq(bm.cardinality(b), 3))
  assert(eq(bm.cardinality(c), 9 - 3))
  assert(eq(bm.hamming(b, c), 9))

end)

test("set/unset many", function ()
  local b = bm.create()
  bm.set(b, 1, 32)
  local c = bm.copy(b)
  bm.unset(c, 17, 32)
  assert(eq(bm.hamming(b, c), 16))
end)

test("set/unset many", function ()
  local b = bm.create()
  assert(eq(bm.tostring(b, 32), "00000000000000000000000000000000"))
end)

test("compress", function ()
  local rand = require("santoku.random")
  local originals = {}
  local n_iterations = 10;
  local n_docs = 1000
  local n_cols_full = 784
  local n_cols_reduced = 512
  for i = 1, n_docs do
    originals[i] = bm.create()
    for j = 1, n_cols_full do
      if rand.num() > 0.9 then
        bm.set(originals[i], j)
      end
    end
  end
  local corpus = bm.matrix(originals, n_cols_full)
  local last = utc.time(true)
  local compress = bm.compressor(
    corpus, n_docs, n_cols_full, n_cols_reduced, n_iterations, nil, function (...)
      local now = utc.time(true)
      print(now - last, ...)
      last = now
    end)
  corpus = compress(corpus, n_docs)
  print(">>", bm.tostring(corpus, n_cols_reduced))
end)
