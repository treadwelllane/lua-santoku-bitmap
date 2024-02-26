local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local vdt = require("santoku.validate")
local eq = vdt.isequal
local bm = require("santoku.bitmap")

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
