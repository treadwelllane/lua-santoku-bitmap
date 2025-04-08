#include "lua.h"
#include "lauxlib.h"
#include "roaring.h"
#include "roaring.c"
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define MT "santoku_bitmap"

static inline void tk_lua_callmod (
  lua_State *L,
  int nargs,
  int nret,
  const char *smod,
  const char *sfn
) {
  lua_getglobal(L, "require"); // arg req
  lua_pushstring(L, smod); // arg req smod
  lua_call(L, 1, 1); // arg mod
  lua_pushstring(L, sfn); // args mod sfn
  lua_gettable(L, -2); // args mod fn
  lua_remove(L, -2); // args fn
  lua_insert(L, - nargs - 1); // fn args
  lua_call(L, nargs, nret); // results
}

static inline int tk_error (
  lua_State *L,
  const char *label,
  int err
) {
  lua_pushstring(L, label);
  lua_pushstring(L, strerror(err));
  tk_lua_callmod(L, 2, 0, "santoku.error", "error");
  return 1;
}

static inline unsigned int tk_lua_checkunsigned (lua_State *L, int i)
{
  lua_Integer l = luaL_checkinteger(L, i);
  if (l < 0)
    luaL_error(L, "value can't be negative");
  if (l > UINT_MAX)
    luaL_error(L, "value is too large");
  return (unsigned int) l;
}

static inline void *tk_malloc (
  lua_State *L,
  size_t s
) {
  void *p = malloc(s);
  if (!p) {
    tk_error(L, "malloc failed", ENOMEM);
    return NULL;
  } else {
    return p;
  }
}

static roaring64_bitmap_t *peek (lua_State *L, int i)
{
  return *((roaring64_bitmap_t **) luaL_checkudata(L, i, MT));
}

static int tk_bitmap_destroy (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  roaring64_bitmap_free(bm);
  return 1;
}

static int tk_bitmap_create (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm = roaring64_bitmap_create();
  if (bm == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bmp = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *)); // t, n, b
  *bmp = bm;
  luaL_getmetatable(L, MT); // t, n, b, mt
  lua_setmetatable(L, -2); // t, n, b
  if (lua_type(L, 1) != LUA_TNIL)
  {
    luaL_checktype(L, 1, LUA_TTABLE);
    lua_Integer n = luaL_checkinteger(L, 2);
    for (lua_Integer i = 1; i <= n; i ++) {
      lua_pushinteger(L, i); // t, n, b, i
      lua_gettable(L, -4); // t, n, b, bit
      if (lua_toboolean(L, -1))
        roaring64_bitmap_add(bm, (uint64_t) i - 1);
      lua_pop(L, 1); // t, n, b
    }
  }
  return 1;
}

static int tk_bitmap_get (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  lua_pushboolean(L, roaring64_bitmap_contains(bm, (uint64_t) bit));
  return 1;
}

static int tk_bitmap_set (lua_State *L)
{
  lua_settop(L, 3);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  if (lua_type(L, 3) != LUA_TNIL) {
    lua_Integer until = luaL_checkinteger(L, 3);
    until --;
    if (until < bit)
      luaL_error(L, "end index must be greater than start index");
    roaring64_bitmap_add_range_closed(bm, (uint64_t) bit, (uint64_t) until);
  } else {
    roaring64_bitmap_add(bm, (uint64_t) bit);
  }
  return 0;
}

static int tk_bitmap_unset (lua_State *L)
{
  lua_settop(L, 3);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  if (lua_type(L, 3) != LUA_TNIL) {
    lua_Integer until = luaL_checkinteger(L, 3);
    until --;
    if (until < bit)
      luaL_error(L, "end index must be greater than start index");
    roaring64_bitmap_remove_range_closed(bm, (uint64_t) bit, (uint64_t) until);
  } else {
    roaring64_bitmap_remove(bm, (uint64_t) bit);
  }
  return 0;
}

static int tk_bitmap_cardinality (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, (lua_Integer) roaring64_bitmap_get_cardinality(bm));
  return 1;
}

static int tk_bitmap_minimum (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, (lua_Integer) roaring64_bitmap_minimum(bm) + 1);
  return 1;
}

static int tk_bitmap_maximum (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, (lua_Integer) roaring64_bitmap_maximum(bm) + 1);
  return 1;
}

static int tk_bitmap_clear (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  roaring64_bitmap_remove_range_closed(bm, 0, UINT64_MAX);
  return 0;
}

typedef struct {
  uint32_t *raw;
} raw_state_t;

static bool tk_bitmap_raw_iter (uint64_t val, void *statepv)
{
  raw_state_t *statep = (raw_state_t *) statepv;
  uint32_t byte = (uint32_t) (val / (sizeof(uint32_t) * CHAR_BIT));
  uint32_t bit = (uint32_t) (val % (sizeof(uint32_t) * CHAR_BIT));
  statep->raw[byte] |= 1 << bit;
  return true;
}

static int tk_bitmap_raw (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm = peek(L, 1);
  uint64_t bits;
  if (lua_type(L, 2) != LUA_TNIL) {
    lua_Integer bs = luaL_checkinteger(L, 2);
    if (bs < 0) {
      return luaL_error(L, "number of bits can't be negative");
    } else if (bits > UINT64_MAX) {
      return luaL_error(L, "number of bits can't be greater than UINT64_MAX");
    } else {
      bits = (uint64_t) bs;
    }
  } else if (roaring64_bitmap_get_cardinality(bm) == 0) {
    bits = 0;
  } else {
    bits = roaring64_bitmap_maximum(bm) + 1;
  }
  uint32_t chunks = bits == 0 ? 1 : (bits - 1) / (sizeof(uint32_t) * CHAR_BIT) + 1;
  raw_state_t state;
  state.raw = tk_malloc(L, sizeof(uint32_t) * chunks);
  memset(state.raw, 0, sizeof(uint32_t) * chunks);
  roaring64_bitmap_iterate(bm, tk_bitmap_raw_iter, &state);
  if (state.raw == NULL)
    luaL_error(L, "error in malloc");
  lua_pushlstring(L, (char *) state.raw, sizeof(uint32_t) * chunks);
  free(state.raw);
  return 1;
}

static int tk_bitmap_from_raw (lua_State *L)
{
  lua_settop(L, 1);
  luaL_checktype(L, 1, LUA_TSTRING);
  size_t size;
  const char *raw = lua_tolstring(L, 1, &size);
  roaring64_bitmap_t *bm = roaring64_bitmap_create();
  if (bm == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bmp = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *)); // s, b
  *bmp = bm;
  luaL_getmetatable(L, MT); // s, b, mt
  lua_setmetatable(L, -2); // s, b
  for (size_t i = 0; i < size; i ++)
    for (size_t j = 0; j < CHAR_BIT; j ++)
      if (raw[i] & (1 << j))
        roaring64_bitmap_add(bm, i * CHAR_BIT + j);
  return 1;
}

static bool tk_bitmap_bits_iter (uint64_t val, void *statepv)
{
  lua_State *L = (lua_State *) statepv; // bm t bits n
  luaL_checktype(L, -3, LUA_TTABLE);
  lua_Integer bits = luaL_checkinteger(L, -2);
  if (bits != -1 && val + 1 > (uint64_t) bits)
    return false;
  lua_Integer n = luaL_checkinteger(L, -1); // bm t bits n
  lua_pushinteger(L, (lua_Integer) val + 1); // bm t bits n v
  lua_settable(L, -4); // bm t bits
  lua_pushinteger(L, n + 1); // bm t bits n
  return val + 1 < (uint64_t) bits;
}

static int tk_bitmap_bits (lua_State *L)
{
  lua_settop(L, 2); // bm bits
  roaring64_bitmap_t *bm = peek(L, 1);
  if (lua_type(L, 2) == LUA_TNIL) {
    lua_pushinteger(L, -1); // bm nil bits
    lua_replace(L, 2); // bm bits
  }
  lua_newtable(L); // bm bits t
  lua_insert(L, 2); // bm t bits
  lua_pushinteger(L, 1); // bm t bits n
  roaring64_bitmap_iterate(bm, tk_bitmap_bits_iter, L); // bm t bits n
  lua_pop(L, 2);
  return 1;
}

static int tk_bitmap_copy (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  roaring64_bitmap_t *bm1 = roaring64_bitmap_copy(bm0);
  if (bm1 == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bmp = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *)); // s, b
  *bmp = bm1;
  luaL_getmetatable(L, MT); // s, b, mt
  lua_setmetatable(L, -2); // s, b
  return 1;
}

static int tk_bitmap_and (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  roaring64_bitmap_t *bm1 = peek(L, 2);
  roaring64_bitmap_and_inplace(bm0, bm1);
  return 0;
}

static int tk_bitmap_equals (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  roaring64_bitmap_t *bm1 = peek(L, 2);
  lua_pushboolean(L, roaring64_bitmap_equals(bm0, bm1));
  return 1;
}

static int tk_bitmap_or (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  roaring64_bitmap_t *bm1 = peek(L, 2);
  roaring64_bitmap_or_inplace(bm0, bm1);
  return 0;
}

static int tk_bitmap_xor (lua_State *L)
{
  lua_settop(L, 2);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  roaring64_bitmap_t *bm1 = peek(L, 2);
  roaring64_bitmap_xor_inplace(bm0, bm1);
  return 0;
}

typedef struct {
  uint64_t n;
  roaring64_bitmap_t *bm;
} extend_state_t;

static bool tk_bitmap_extend_iter (uint64_t val, void *statepv)
{
  extend_state_t *statep = (extend_state_t *) statepv;
  roaring64_bitmap_add(statep->bm, val + statep->n);
  return true;
}

static int tk_bitmap_extend (lua_State *L)
{
  lua_settop(L, 3);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  roaring64_bitmap_t *bm1 = peek(L, 2);
  lua_Integer n = luaL_checkinteger(L, 3);
  n --;
  if (n < 0)
    luaL_error(L, "extension starting index must be greater than 0");
  extend_state_t state;
  state.n = (uint64_t) n;
  state.bm = bm0;
  roaring64_bitmap_iterate(bm1, tk_bitmap_extend_iter, &state);
  return 0;
}

static int tk_bitmap_flip (lua_State *L)
{
  lua_settop(L, 3);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  if (lua_type(L, 3) != LUA_TNIL) {
    lua_Integer until = luaL_checkinteger(L, 3);
    until --;
    if (until < bit)
      luaL_error(L, "end index must be greater than start index");
    roaring64_bitmap_flip_inplace(bm, (uint64_t) bit, (uint64_t) until + 1);
  } else {
    roaring64_bitmap_flip_inplace(bm, 0, (uint64_t) bit + 1);
  }
  return 0;
}

static int tk_bitmap_flip_interleave (lua_State *L)
{
  lua_settop(L, 3);
  roaring64_bitmap_t *bm0 = peek(L, 1);
  uint64_t n = tk_lua_checkunsigned(L, 2);
  uint64_t step = tk_lua_checkunsigned(L, 3);
  roaring64_bitmap_t *bm1 = roaring64_bitmap_create();
  if (bm1 == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bm1p = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *)); // s, b
  *bm1p = bm1;
  luaL_getmetatable(L, MT); // s, b, mt
  lua_setmetatable(L, -2); // s, b
  // TODO: is there a faster way to do this with set operations?
  for (uint64_t i = 0; i < n; i ++)
    for (uint64_t j = 0; j < step; j ++)
      if (roaring64_bitmap_contains(bm0, i * step + j))
        roaring64_bitmap_add(bm1,
          i * 2 * step + 0 * step + j);
      else
        roaring64_bitmap_add(bm1,
          i * 2 * step + 1 * step + j);
  return 1;
}

static luaL_Reg fns[] =
{
  { "create", tk_bitmap_create },
  { "copy", tk_bitmap_copy },
  { "destroy", tk_bitmap_destroy },
  { "set", tk_bitmap_set },
  { "get", tk_bitmap_get },
  { "unset", tk_bitmap_unset },
  { "cardinality", tk_bitmap_cardinality },
  { "minimum", tk_bitmap_minimum },
  { "maximum", tk_bitmap_maximum },
  { "clear", tk_bitmap_clear },
  { "bits", tk_bitmap_bits },
  { "raw", tk_bitmap_raw },
  { "from_raw", tk_bitmap_from_raw },
  { "equals", tk_bitmap_equals },
  { "and", tk_bitmap_and },
  { "or", tk_bitmap_or },
  { "xor", tk_bitmap_xor },
  { "flip", tk_bitmap_flip },
  { "flip_interleave", tk_bitmap_flip_interleave },
  { "extend", tk_bitmap_extend },
  { "bits", tk_bitmap_bits },
  { NULL, NULL }
};

int luaopen_santoku_bitmap_capi (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, fns); // t
  lua_pushinteger(L, sizeof(uint32_t) * CHAR_BIT); // t i
  lua_setfield(L, -2, "chunk_bits"); // t
  luaL_newmetatable(L, MT); // t mt
  lua_pushcfunction(L, tk_bitmap_destroy); // t mt fn
  lua_setfield(L, -2, "__gc"); // t mt
  lua_pop(L, 1); // t
  return 1;
}
