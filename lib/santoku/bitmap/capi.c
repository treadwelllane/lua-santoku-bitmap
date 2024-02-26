#include "lua.h"
#include "lauxlib.h"
#include "roaring.c"
#include <string.h>

#define TK_BITMAP_MT "santoku_bitmap"

roaring_bitmap_t *tk_bitmap_peek (lua_State *L, int i)
{
  return *((roaring_bitmap_t **) luaL_checkudata(L, i, TK_BITMAP_MT));
}

int tk_bitmap_destroy (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
  roaring_bitmap_free(bm);
  return 1;
}

int tk_bitmap_create (lua_State *L)
{
  lua_settop(L, 0);
  roaring_bitmap_t *bm = roaring_bitmap_create();
  if (bm == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring_bitmap_t **bmp = (roaring_bitmap_t **) lua_newuserdata(L, sizeof(roaring_bitmap_t *));
  *bmp = bm;
  luaL_getmetatable(L, TK_BITMAP_MT);
  lua_setmetatable(L, -2);
  return 1;
}

int tk_bitmap_get (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  lua_pushboolean(L, roaring_bitmap_contains(bm, bit));
  return 1;
}

int tk_bitmap_set (lua_State *L)
{
  lua_settop(L, 3);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
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

int tk_bitmap_unset (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
  lua_Integer bit = luaL_checkinteger(L, 2);
  bit --;
  if (bit < 0)
    luaL_error(L, "bit index must be greater than zero");
  roaring_bitmap_remove(bm, bit);
  return 0;
}

int tk_bitmap_cardinality (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
  lua_pushinteger(L, roaring_bitmap_get_cardinality(bm));
  return 1;
}

int tk_bitmap_clear (lua_State *L)
{
  lua_settop(L, 1);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
  roaring_bitmap_clear(bm);
  return 0;
}

typedef struct {
  unsigned int *raw;
} tk_bitmap_raw_state_t;

bool tk_bitmap_raw_iter (uint32_t val, void *statepv)
{
  tk_bitmap_raw_state_t *statep = (tk_bitmap_raw_state_t *) statepv;
  unsigned int byte = val / (sizeof(unsigned int) * CHAR_BIT);
  unsigned int bit = val % (sizeof(unsigned int) * CHAR_BIT);
  statep->raw[byte] |= 1 << bit;
  return true;
}

int tk_bitmap_raw (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm = tk_bitmap_peek(L, 1);
  size_t chunks;
  if (lua_type(L, 2) != LUA_TNIL) {
    lua_Integer bits = luaL_checkinteger(L, 2);
    if (bits == 0) {
      chunks = 0;
    } else if (bits < 0) {
      return luaL_error(L, "number of bits can't be negative");
    } else if (bits > UINT32_MAX) {
      return luaL_error(L, "number of bits can't be greater than UINT32_MAX");
    } else {
      chunks = bits / (sizeof(unsigned int) * CHAR_BIT) + 1;
    }
  } else if (roaring_bitmap_get_cardinality(bm) == 0) {
    chunks = 0;
  } else {
    uint32_t max = roaring_bitmap_maximum(bm);
    chunks = max / (sizeof(unsigned int) * CHAR_BIT) + 1;
  }
  tk_bitmap_raw_state_t state;
  state.raw = malloc(sizeof(unsigned int) * chunks);
  memset(state.raw, 0, sizeof(unsigned int) * chunks);
  roaring_iterate(bm, tk_bitmap_raw_iter, &state);
  if (state.raw == NULL)
    luaL_error(L, "error in malloc");
  lua_pushlstring(L, (char *) state.raw, sizeof(unsigned int) * chunks);
  free(state.raw);
  return 1;
}

int tk_bitmap_tostring (lua_State *L)
{
  lua_settop(L, 2);
  tk_bitmap_raw(L);
  size_t size_c;
  const char *raw_c = luaL_checklstring(L, -1, &size_c);
  size_t size_u = size_c ? size_c / sizeof(unsigned int) : 0;
  unsigned int *raw_u = (unsigned int *) raw_c;
  luaL_Buffer buf;
  luaL_buffinit(L, &buf);
  for (size_t i = 0; i < size_u; i ++)
    for (unsigned int c = 0; c < sizeof(unsigned int) * CHAR_BIT; c ++)
      luaL_addchar(&buf, (raw_u[i] & (1 << c)) ? '1' : '0');
  luaL_pushresult(&buf);
  return 1;
}

int tk_bitmap_and (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = tk_bitmap_peek(L, 1);
  roaring_bitmap_t *bm1 = tk_bitmap_peek(L, 2);
  roaring_bitmap_and_inplace(bm0, bm1);
  return 0;
}

int tk_bitmap_equals (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = tk_bitmap_peek(L, 1);
  roaring_bitmap_t *bm1 = tk_bitmap_peek(L, 2);
  lua_pushboolean(L, roaring_bitmap_equals(bm0, bm1));
  return 1;
}

int tk_bitmap_or (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = tk_bitmap_peek(L, 1);
  roaring_bitmap_t *bm1 = tk_bitmap_peek(L, 2);
  roaring_bitmap_or_inplace(bm0, bm1);
  return 0;
}

int tk_bitmap_xor (lua_State *L)
{
  lua_settop(L, 2);
  roaring_bitmap_t *bm0 = tk_bitmap_peek(L, 1);
  roaring_bitmap_t *bm1 = tk_bitmap_peek(L, 2);
  roaring_bitmap_xor_inplace(bm0, bm1);
  return 0;
}

typedef struct {
  uint32_t n;
  roaring_bitmap_t *bm;
} tk_bitmap_extend_state_t;

bool tk_bitmap_extend_iter (uint32_t val, void *statepv)
{
  tk_bitmap_extend_state_t *statep = (tk_bitmap_extend_state_t *) statepv;
  roaring_bitmap_add(statep->bm, val + statep->n);
  return true;
}

int tk_bitmap_extend (lua_State *L)
{
  lua_settop(L, 3);
  roaring_bitmap_t *bm0 = tk_bitmap_peek(L, 1);
  roaring_bitmap_t *bm1 = tk_bitmap_peek(L, 2);
  lua_Integer n = luaL_checkinteger(L, 3);
  n --;
  if (n < 0)
    luaL_error(L, "extension starting index must be greater than 0");
  tk_bitmap_extend_state_t state;
  state.n = n;
  state.bm = bm0;
  roaring_iterate(bm1, tk_bitmap_extend_iter, &state);
  return 0;
}

luaL_Reg tk_bitmap_fns[] =
{
  { "create", tk_bitmap_create },
  { "destroy", tk_bitmap_destroy },
  { "set", tk_bitmap_set },
  { "get", tk_bitmap_get },
  { "unset", tk_bitmap_unset },
  { "cardinality", tk_bitmap_cardinality },
  { "clear", tk_bitmap_clear },
  { "raw", tk_bitmap_raw },
  { "tostring", tk_bitmap_tostring },
  { "equals", tk_bitmap_equals },
  { "and", tk_bitmap_and },
  { "or", tk_bitmap_or },
  { "xor", tk_bitmap_xor },
  { "extend", tk_bitmap_extend },
  { NULL, NULL }
};

int luaopen_santoku_bitmap_capi (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_bitmap_fns); // t
  luaL_newmetatable(L, TK_BITMAP_MT); // t mt
  lua_pushcfunction(L, tk_bitmap_destroy); // t mt fn
  lua_setfield(L, -2, "__gc"); // t mt
  lua_pop(L, 1); // t
  return 1;
}
