#include "lua.h"
#include "lauxlib.h"
#include "roaring.c"
#include <string.h>

#define MT "santoku_bitmap"

static roaring_bitmap_t *peek (lua_State *L, int i)
{
  return *((roaring_bitmap_t **) luaL_checkudata(L, i, MT));
}

static int destroy (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = peek(L, 1);
  roaring_bitmap_free(bm);
  return 1;
}

static int create (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = roaring_bitmap_create();
  if (bm == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring_bitmap_t **bmp = (roaring_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring_bitmap_t *)); // t, n, b
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
        roaring_bitmap_add(bm, i - 1);
      lua_pop(L, 1); // t, n, b
    }
  }
  return 1;
}

static int get (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  lua_pushboolean(L, roaring_bitmap_contains(bm, bit));
  return 1;
}

static int set (lua_State *L)
{
  lua_settop(L, 3);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  if (lua_type(L, 3) != LUA_TNIL) {
    lua_Integer until = luaL_checkinteger(L, 3);
    until --;
    if (until < bit)
      luaL_error(L, "end index must be greater than start index");
    roaring_bitmap_add_range_closed(bm, bit, until);
  } else {
    roaring_bitmap_add(bm, bit);
  }
  return 0;
}

static int unset (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  roaring_bitmap_remove(bm, bit);
  return 0;
}

static int cardinality (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, roaring_bitmap_get_cardinality(bm));
  return 1;
}

static int minimum (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, roaring_bitmap_minimum(bm) + 1);
  return 1;
}

static int maximum (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, roaring_bitmap_maximum(bm) + 1);
  return 1;
}

static int clear (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = peek(L, 1);
  roaring_bitmap_clear(bm);
  return 0;
}

typedef struct {
  uint32_t *raw;
} raw_state_t;

static bool raw_iter (uint32_t val, void *statepv)
{
  raw_state_t *statep = (raw_state_t *) statepv;
  uint32_t byte = val / (sizeof(uint32_t) * CHAR_BIT);
  uint32_t bit = val % (sizeof(uint32_t) * CHAR_BIT);
  statep->raw[byte] |= 1 << bit;
  return true;
}

static int raw (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_Integer bits;
  if (lua_type(L, 2) != LUA_TNIL) {
    bits = luaL_checkinteger(L, 2);
    if (bits < 0) {
      return luaL_error(L, "number of bits can't be negative");
    } else if (bits > UINT32_MAX) {
      return luaL_error(L, "number of bits can't be greater than UINT32_MAX");
    }
  } else if (roaring_bitmap_get_cardinality(bm) == 0) {
    bits = 0;
  } else {
    bits = roaring_bitmap_maximum(bm) + 1;
  }
  uint32_t chunks = (bits - 1) / (sizeof(uint32_t) * CHAR_BIT) + 1;
  raw_state_t state;
  state.raw = malloc(sizeof(uint32_t) * chunks);
  memset(state.raw, 0, sizeof(uint32_t) * chunks);
  roaring_iterate(bm, raw_iter, &state);
  if (state.raw == NULL)
    luaL_error(L, "error in malloc");
  lua_pushlstring(L, (char *) state.raw, sizeof(uint32_t) * chunks);
  free(state.raw);
  return 1;
}

static int from_raw (lua_State *L)
{
  lua_settop(L, 1);
  luaL_checktype(L, 1, LUA_TSTRING);
  size_t size;
  const char *raw = lua_tolstring(L, 1, &size);
  roaring_bitmap_t *bm = roaring_bitmap_create();
  if (bm == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring_bitmap_t **bmp = (roaring_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring_bitmap_t *)); // s, b
  *bmp = bm;
  luaL_getmetatable(L, MT); // s, b, mt
  lua_setmetatable(L, -2); // s, b
  for (size_t i = 0; i < size; i ++)
    for (size_t j = 0; j < CHAR_BIT; j ++)
      if (raw[i] & (1 << j))
        roaring_bitmap_add(bm, i * CHAR_BIT + j);
  return 1;
}

static int copy (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm0 = peek(L, 1);
  roaring_bitmap_t *bm1 = roaring_bitmap_copy(bm0);
  roaring_bitmap_t **bmp = (roaring_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring_bitmap_t *)); // s, b
  *bmp = bm1;
  luaL_getmetatable(L, MT); // s, b, mt
  lua_setmetatable(L, -2); // s, b
  return 1;
}

static int tostring (lua_State *L)
{
  lua_settop(L, 2);
  raw(L);
  size_t size_c;
  const char *raw_c = luaL_checklstring(L, -1, &size_c);
  size_t size_u = size_c ? size_c / sizeof(uint32_t) : 0;
  uint32_t *raw_u = (uint32_t *) raw_c;
  luaL_Buffer buf;
  luaL_buffinit(L, &buf);
  for (size_t i = 0; i < size_u; i ++)
    for (uint32_t c = 0; c < sizeof(uint32_t) * CHAR_BIT; c ++)
      luaL_addchar(&buf, (raw_u[i] & (1 << c)) ? '1' : '0');
  luaL_pushresult(&buf);
  return 1;
}

static int and (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = peek(L, 1);
  roaring_bitmap_t *bm1 = peek(L, 2);
  roaring_bitmap_and_inplace(bm0, bm1);
  return 0;
}

static int equals (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = peek(L, 1);
  roaring_bitmap_t *bm1 = peek(L, 2);
  lua_pushboolean(L, roaring_bitmap_equals(bm0, bm1));
  return 1;
}

static int or (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = peek(L, 1);
  roaring_bitmap_t *bm1 = peek(L, 2);
  roaring_bitmap_or_inplace(bm0, bm1);
  return 0;
}

static int xor (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = peek(L, 1);
  roaring_bitmap_t *bm1 = peek(L, 2);
  roaring_bitmap_xor_inplace(bm0, bm1);
  return 0;
}

typedef struct {
  uint32_t n;
  roaring_bitmap_t *bm;
} extend_state_t;

static bool extend_iter (uint32_t val, void *statepv)
{
  extend_state_t *statep = (extend_state_t *) statepv;
  roaring_bitmap_add(statep->bm, val + statep->n);
  return true;
}

static int extend (lua_State *L)
{
  lua_settop(L, 3);
  roaring_bitmap_t *bm0 = peek(L, 1);
  roaring_bitmap_t *bm1 = peek(L, 2);
  lua_Integer n = luaL_checkinteger(L, 3);
  n --;
  if (n < 0)
    luaL_error(L, "extension starting index must be greater than 0");
  extend_state_t state;
  state.n = n;
  state.bm = bm0;
  roaring_iterate(bm1, extend_iter, &state);
  return 0;
}

static int flip (lua_State *L)
{
  lua_settop(L, 3);
  roaring_bitmap_t *bm = peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  if (lua_type(L, 3) != LUA_TNIL) {
    lua_Integer until = luaL_checkinteger(L, 3);
    until --;
    if (until < bit)
      luaL_error(L, "end index must be greater than start index");
    roaring_bitmap_flip(bm, bit, until);
  } else {
    roaring_bitmap_flip(bm, 0, bit);
  }
  return 0;
}

static luaL_Reg fns[] =
{
  { "create", create },
  { "copy", copy },
  { "destroy", destroy },
  { "set", set },
  { "get", get },
  { "unset", unset },
  { "cardinality", cardinality },
  { "minimum", minimum },
  { "maximum", maximum },
  { "clear", clear },
  { "raw", raw },
  { "from_raw", from_raw },
  { "tostring", tostring },
  { "equals", equals },
  { "and", and },
  { "or", or },
  { "xor", xor },
  { "flip", flip },
  { "extend", extend },
  { NULL, NULL }
};

int luaopen_santoku_bitmap_capi (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, fns); // t
  lua_pushinteger(L, sizeof(uint32_t) * CHAR_BIT); // t i
  lua_setfield(L, -2, "chunk_bits"); // t
  luaL_newmetatable(L, MT); // t mt
  lua_pushcfunction(L, destroy); // t mt fn
  lua_setfield(L, -2, "__gc"); // t mt
  lua_pop(L, 1); // t
  return 1;
}
