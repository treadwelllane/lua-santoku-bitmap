local test = require("santoku.test")
local utc = require("santoku.utc")
local str = require("santoku.string")
local err = require("santoku.error")
local fs = require("santoku.fs")
local assert = err.assert
local vdt = require("santoku.validate")
local eq = vdt.isequal
local bm = require("santoku.bitmap")
local mtx = require("santoku.matrix")
local bmc = require("santoku.bitmap.compressor")

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

test("flip_interleave", function ()
  local b = bm.create()
  bm.set(b, 0 + 1)
  bm.set(b, 0 + 4)
  bm.set(b, 0 + 5)
  bm.set(b, 0 + 6)
  bm.set(b, 6 + 3)
  bm.set(b, 6 + 4)
  b = bm.flip_interleave(b, 2, 6)
  assert(eq(bm.tostring(b, 12), "10011101100000110011001100000000"))
end)

test("compress", function ()
  local rand = require("santoku.random")
  local originals = {}
  local labels = {}
  local n_iterations = 20
  local n_docs = 10000
  local n_cols_full = 1786
  local n_cols_reduced = 128
  local n_threads = nil
  for i = 1, n_docs do
    originals[i] = bm.create()
    labels[i] = rand.num() > 0.5 and 1 or 0
    for j = 1, n_cols_full do
      if rand.num() > 0.95 then
        bm.set(originals[i], j)
      end
    end
  end
  local corpus_original = bm.matrix(originals, n_cols_full)
  labels = mtx.raw(mtx.create(labels), nil, nil, "u32")
  print()
  print(require("santoku.serialize")(bm.top_mi(
    corpus_original, labels, n_docs, n_cols_full, 2, 100), true))
  print(require("santoku.serialize")(bm.top_chi2(
    corpus_original, labels, n_docs, n_cols_full, 2, 100), true))
  local compressor = bmc.create({
    visible = n_cols_full,
    hidden = n_cols_reduced,
    threads = n_threads,
  })
  local stopwatch = utc.stopwatch(0.1)
  compressor.train({
    corpus = corpus_original,
    samples = n_docs,
    iterations = n_iterations,
    each = function (i, c, d)
      local duration = stopwatch()
      str.printf(" %3d   %3.4f   %6.4f   %6.4f\n", i, duration, c, d)
      return true
    end
  })
  local corpus_compressed = compressor.compress(corpus_original, n_docs)
  local corpus_compressed0 = compressor.compress(corpus_original, n_docs)
  assert(bm.equals(corpus_compressed, corpus_compressed))
  assert(bm.equals(corpus_compressed, corpus_compressed0))
  fs.rm(".tmp.compressor.bin", true)
  compressor.persist(".tmp.compressor.bin")
  compressor.destroy()
  compressor = bmc.load(".tmp.compressor.bin", 3)
  local corpus_compressed1 = compressor.compress(corpus_original, n_docs)
  assert(bm.equals(corpus_compressed, corpus_compressed))
  assert(bm.equals(corpus_compressed, corpus_compressed1))
end)
