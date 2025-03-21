#include "lua.h"
#include "lauxlib.h"
#include "roaring.h"
#include "roaring.c"
#include <string.h>
#include <math.h>
#include <float.h>

#define MT "santoku_bitmap"

static inline int tk_lua_absindex (lua_State *L, int i)
{
  if (i < 0 && i > LUA_REGISTRYINDEX)
    i += lua_gettop(L) + 1;
  return i;
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

static inline unsigned int tk_lua_optunsigned (lua_State *L, int i, unsigned int def)
{
  if (lua_type(L, i) < 1)
    return def;
  return tk_lua_checkunsigned(L, i);
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
        roaring64_bitmap_add(bm, i - 1);
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
  lua_pushboolean(L, roaring64_bitmap_contains(bm, bit));
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
    roaring64_bitmap_add_range_closed(bm, bit, until);
  } else {
    roaring64_bitmap_add(bm, bit);
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
    roaring64_bitmap_remove_range_closed(bm, bit, until);
  } else {
    roaring64_bitmap_remove(bm, bit);
  }
  return 0;
}

static int tk_bitmap_cardinality (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, roaring64_bitmap_get_cardinality(bm));
  return 1;
}

static int tk_bitmap_minimum (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, roaring64_bitmap_minimum(bm) + 1);
  return 1;
}

static int tk_bitmap_maximum (lua_State *L)
{
  lua_settop(L, 1);
  roaring64_bitmap_t *bm = peek(L, 1);
  lua_pushinteger(L, roaring64_bitmap_maximum(bm) + 1);
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
  state.raw = malloc(sizeof(uint32_t) * chunks);
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
  if (val + 1 > bits)
    return false;
  lua_Integer n = luaL_checkinteger(L, -1); // bm t bits n
  lua_pushinteger(L, val + 1); // bm t bits n v
  lua_settable(L, -4); // bm t bits
  lua_pushinteger(L, n + 1); // bm t bits n
  return val + 1 < bits;
}

static int tk_bitmap_bits (lua_State *L)
{
  lua_settop(L, 2); // bm bits
  roaring64_bitmap_t *bm = peek(L, 1);
  if (lua_type(L, 2) == LUA_TNIL) {
    lua_pushinteger(L, UINT64_MAX); // bm nil bits
    lua_replace(L, 2); // bm bits
  }
  lua_newtable(L); // bm bits t
  lua_insert(L, 2); // bm t bits
  lua_pushinteger(L, 1); // bm t bits n
  roaring64_bitmap_iterate(bm, tk_bitmap_bits_iter, L); // bm t bits n
  lua_pushnil(L); // bm t bits n nil
  lua_settable(L, -4); // bm t bits
  lua_pop(L, 1);
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
  state.n = n;
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
    roaring64_bitmap_flip_inplace(bm, bit, until + 1);
  } else {
    roaring64_bitmap_flip_inplace(bm, 0, bit + 1);
  }
  return 0;
}

#define MT_COMPRESSOR "santoku_bitmap_compressor"

typedef struct {
  double *alpha;
  double *log_marg;
  double *log_p_y;
  double *log_p_y_given_x_unnorm;
  double *maxmis;
  double *mis;
  double *sums;
  double *sumsa;
  double *p_x;
  double *entropy_x;
  double *p_y_given_x;
  double *pseudo_counts;
  double *smis;
  double *tc_history;
  double *tcs;
  double eps;
  double lam;
  double tmin;
  double ttc;
  unsigned int n_visible;
  unsigned int n_hidden;
} tk_compressor_t;

static tk_compressor_t *peek_compressor (lua_State *L, int i)
{
  return (tk_compressor_t *) luaL_checkudata(L, i, MT_COMPRESSOR);
}

static inline void tk_compressor_shrink (tk_compressor_t *compressor)
{
  if (compressor->entropy_x) {
    free(compressor->entropy_x);
    compressor->entropy_x = NULL;
  }
  if (compressor->maxmis) {
    free(compressor->maxmis);
    compressor->maxmis = NULL;
  }
  if (compressor->mis) {
    free(compressor->mis);
    compressor->mis = NULL;
  }
  if (compressor->sums) {
    free(compressor->sums);
    compressor->sums = NULL;
  };
  if (compressor->p_x) {
    free(compressor->p_x);
    compressor->p_x = NULL;
  }
  if (compressor->pseudo_counts) {
    free(compressor->pseudo_counts);
    compressor->pseudo_counts = NULL;
  }
  if (compressor->smis) {
    free(compressor->smis);
    compressor->smis = NULL;
  }
  if (compressor->tc_history) {
    free(compressor->tc_history);
    compressor->tc_history = NULL;
  }
  if (compressor->tcs) {
    free(compressor->tcs);
    compressor->tcs = NULL;
  }
}

static int tk_compressor_destroy (lua_State *L)
{
  lua_settop(L, 1);
  tk_compressor_t *compressor = peek_compressor(L, 1);
  tk_compressor_shrink(compressor);
  if (compressor->alpha) {
    free(compressor->alpha);
    compressor->alpha = NULL;
  }
  if (compressor->log_marg) {
    free(compressor->log_marg);
    compressor->log_marg = NULL;
  }
  if (compressor->log_p_y) {
    free(compressor->log_p_y);
    compressor->log_p_y = NULL;
  }
  if (compressor->log_p_y_given_x_unnorm) {
    free(compressor->log_p_y_given_x_unnorm);
    compressor->log_p_y_given_x_unnorm = NULL;
  }
  if (compressor->p_y_given_x) {
    free(compressor->p_y_given_x);
    compressor->p_y_given_x = NULL;
  }
  if (compressor->sumsa) {
    free(compressor->sumsa);
    compressor->sumsa = NULL;
  };
  return 1;
}

static inline int tk_bitmap_compress (lua_State *);

static inline void tk_compressor_logsumexp (
  double *restrict log_p_y_given_x_unnorm,
  double *restrict log_z, // mis
  unsigned int n_hidden,
  unsigned int n_samples
) {
  double *restrict offset0 = log_p_y_given_x_unnorm + 0 * n_hidden * n_samples;
  double *restrict offset1 = log_p_y_given_x_unnorm + 1 * n_hidden * n_samples;
  for (unsigned int i = 0; i < n_hidden * n_samples; i ++) {
    double a = offset0[i];
    double b = offset1[i];
    double max_ab  = (a > b) ? a : b;
    double sum_exp = exp(a - max_ab) + exp(b - max_ab);
    log_z[i] = max_ab + log(sum_exp);
  }
}

static void tk_compressor_normalize_latent (
  double *restrict log_z, // mis
  double *restrict p_y_given_x,
  double *restrict log_p_y_given_x_unnorm,
  unsigned int n_hidden,
  unsigned int n_samples
) {
  tk_compressor_logsumexp(log_p_y_given_x_unnorm, log_z, n_hidden, n_samples);
  double *restrict pyx_offset0 = p_y_given_x + 0 * n_hidden * n_samples;
  double *restrict pyx_offset1 = p_y_given_x + 1 * n_hidden * n_samples;
  double *restrict pyxu_offset0 = log_p_y_given_x_unnorm + 0 * n_hidden * n_samples;
  double *restrict pyxu_offset1 = log_p_y_given_x_unnorm + 1 * n_hidden * n_samples;
  for (unsigned int i = 0; i < n_hidden * n_samples; i ++)
    pyx_offset0[i] = exp(pyxu_offset0[i] - log_z[i]);
  for (unsigned int i = 0; i < n_hidden * n_samples; i ++)
    pyx_offset1[i] = exp(pyxu_offset1[i] - log_z[i]);
}

static inline void tk_compressor_data_stats (
  uint64_t cardinality,
  uint64_t *restrict bits,
  double *restrict p_x,
  double *restrict entropy_x,
  unsigned int n_samples,
  unsigned int n_visible
) {
  for (unsigned int v = 0; v < n_visible; v ++)
    p_x[v] = 0;
  for (uint64_t c = 0; c < cardinality; c ++) // not vectorized, aliasing
    p_x[bits[c] % n_visible] ++;
  for (unsigned int v = 0; v < n_visible; v ++)
    p_x[v] /= (double) n_samples;
  for (unsigned int v = 0; v < n_visible; v ++) {
    double entropy = 0;
    entropy -= p_x[v] * log(p_x[v]);
    entropy -= (1 - p_x[v]) * log(p_x[v]);
    entropy_x[v] = entropy > 0 ? entropy : 1e-10;
  }
}

static inline void tk_compressor_calculate_p_y (
  double *restrict log_p_y,
  double *restrict p_y_given_x,
  unsigned int n_samples,
  unsigned int n_hidden
) {

  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, unsupported outer form
    double *restrict offset0 = p_y_given_x + 0 * n_hidden * n_samples + h * n_samples;
    double *restrict offset1 = p_y_given_x + 1 * n_hidden * n_samples + h * n_samples;
    double pseudo_counts_0 = 0.001;
    double pseudo_counts_1 = 0.001;
    for (unsigned int s = 0; s < n_samples; s ++) {
      pseudo_counts_0 += offset0[s];
      pseudo_counts_1 += offset1[s];
    }
    double sum_pseudo_counts = pseudo_counts_0 + pseudo_counts_1;
    log_p_y[h * 2 + 0] = log(pseudo_counts_0) - log(sum_pseudo_counts);
    log_p_y[h * 2 + 1] = log(pseudo_counts_1) - log(sum_pseudo_counts);
  }
}

static void tk_compressor_calculate_p_y_xi (
  uint64_t cardinality,
  uint64_t *restrict bits,
  double *restrict log_marg,
  double *restrict pseudo_counts,
  double *restrict log_p_y,
  double *restrict p_y_given_x,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  for (unsigned int i = 0; i < 2 * 2 * n_hidden * n_visible; i ++)
    pseudo_counts[i] = 0;
  double *restrict pc00 = pseudo_counts + 0 * n_hidden * n_visible;
  double *restrict pc01 = pseudo_counts + 1 * n_hidden * n_visible;
  double *restrict pc10 = pseudo_counts + 2 * n_hidden * n_visible;
  double *restrict pc11 = pseudo_counts + 3 * n_hidden * n_visible;
  double *restrict lm00 = log_marg + 0 * n_hidden * n_visible;
  double *restrict lm01 = log_marg + 1 * n_hidden * n_visible;
  double *restrict lm10 = log_marg + 2 * n_hidden * n_visible;
  double *restrict lm11 = log_marg + 3 * n_hidden * n_visible;
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, 2+ inner loops
    double sum_py0 = 0.0;
    double sum_py1 = 0.0;
    double *restrict hpy0s = p_y_given_x + 0 * n_hidden * n_samples + h * n_samples;
    double *restrict hpy1s = p_y_given_x + 1 * n_hidden * n_samples + h * n_samples;
    for (unsigned int s = 0; s < n_samples; s ++)
      sum_py0 += hpy0s[s];
    for (unsigned int s = 0; s < n_samples; s ++)
      sum_py1 += hpy1s[s];
    double *restrict pc00a = pc00 + h * n_visible;
    double *restrict pc01a = pc01 + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      pc00a[v] += sum_py0;
      pc01a[v] += sum_py1;
    }
  }
  for (unsigned int h = 0; h < n_hidden; h ++) {
    for (unsigned int v = 0; v < n_visible; v ++) {
      lm00[v] = 0;
      lm01[v] = 0;
      lm10[v] = 0;
      lm11[v] = 0;
    }
    for (uint64_t c = 0; c < cardinality; c ++) { // not vectorzed, gather scatter
      uint64_t b = bits[c];
      uint64_t s = b / n_visible;
      uint64_t v = b % n_visible;
      double py0 = p_y_given_x[0 * n_hidden * n_samples + h * n_samples + s];
      double py1 = p_y_given_x[1 * n_hidden * n_samples + h * n_samples + s];
      lm10[v] += py0;
      lm11[v] += py1;
      lm00[v] -= py0;
      lm01[v] -= py1;
    }
    double *restrict pc00a = pc00 + h * n_visible;
    double *restrict pc01a = pc01 + h * n_visible;
    double *restrict pc10a = pc10 + h * n_visible;
    double *restrict pc11a = pc11 + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      pc10a[v] += lm10[v];
      pc11a[v] += lm11[v];
      pc00a[v] += lm00[v];
      pc01a[v] += lm01[v];
    }
  }
  for (unsigned int i = 0; i < 2 * 2 * n_hidden * n_visible; i ++) // splitting at boundaries
    pseudo_counts[i] += 0.001;
  for (unsigned int i = 0; i < n_hidden * n_visible; i ++) { // splitting at boundaries
    double log_total0 = log(pc00[i] + pc01[i]);
    lm00[i] = log(pc00[i]) - log_total0;
    lm01[i] = log(pc01[i]) - log_total0;
  }
  for (unsigned int i = 0; i < n_hidden * n_visible; i ++) { // splitting at boundaries
    double log_total1 = log(pc10[i] + pc11[i]);
    lm10[i] = log(pc10[i]) - log_total1;
    lm11[i] = log(pc11[i]) - log_total1;
  }
  for (unsigned int h = 0; h < n_hidden; h ++) { // not affine
    double *restrict lm00a = lm00 + h * n_visible;
    double *restrict lm01a = lm01 + h * n_visible;
    double *restrict lm10a = lm10 + h * n_visible;
    double *restrict lm11a = lm11 + h * n_visible;
    double lpy0 = log_p_y[h * 2 + 0];
    double lpy1 = log_p_y[h * 2 + 1];
    for (unsigned int v = 0; v < n_visible; v ++) {
      lm00a[v] -= lpy0;
      lm01a[v] -= lpy1;
      lm10a[v] -= lpy0;
      lm11a[v] -= lpy1;
    }
  }
}

static inline void tk_compressor_calculate_mis (
  double *restrict mis,
  double *restrict smis,
  double *restrict pseudo_counts,
  double *restrict log_p_y,
  double *restrict log_marg,
  double *restrict p_x,
  double *restrict entropy_x,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  double *restrict lm00 = log_marg + 0 * n_hidden * n_visible;
  double *restrict lm01 = log_marg + 1 * n_hidden * n_visible;
  double *restrict lm10 = log_marg + 2 * n_hidden * n_visible;
  double *restrict lm11 = log_marg + 3 * n_hidden * n_visible;
  double *restrict vc00 = pseudo_counts + 0 * n_hidden * n_visible;
  double *restrict vc01 = pseudo_counts + 1 * n_hidden * n_visible;
  double *restrict vc10 = pseudo_counts + 2 * n_hidden * n_visible;
  double *restrict vc11 = pseudo_counts + 3 * n_hidden * n_visible;
  double *restrict smis0 = smis + 0 * n_hidden * n_visible;
  double *restrict smis1 = smis + 1 * n_hidden * n_visible;
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, non-affine base
    double *restrict lm00a = lm00 + h * n_visible;
    double *restrict lm01a = lm01 + h * n_visible;
    double *restrict lm10a = lm10 + h * n_visible;
    double *restrict lm11a = lm11 + h * n_visible;
    double *restrict vc00a = vc00 + h * n_visible;
    double *restrict vc01a = vc01 + h * n_visible;
    double *restrict vc10a = vc10 + h * n_visible;
    double *restrict vc11a = vc11 + h * n_visible;
    double lpy0 = log_p_y[h * 2 + 0];
    double lpy1 = log_p_y[h * 2 + 1];
    for (unsigned int v = 0; v < n_visible; v ++) {
      vc00a[v] = exp(lm00a[v] + lpy0);
      vc01a[v] = exp(lm01a[v] + lpy1);
      vc10a[v] = exp(lm10a[v] + lpy0);
      vc11a[v] = exp(lm11a[v] + lpy1);
    }
  }
  for (unsigned int i = 0; i < n_hidden * n_visible; i ++) // not vectorized, splitting dominance boundary
    smis0[i] = vc00[i] * lm00[i] + vc01[i] * lm01[i];
  for (unsigned int i = 0; i < n_hidden * n_visible; i ++)
    smis1[i] = vc10[i] * lm10[i] + vc11[i] * lm11[i];
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, dominance boundary
    double *restrict smis0a = smis0 + h * n_visible;
    double *restrict smis1a = smis1 + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      double weighted_sum = smis0a[v] * (1 - p_x[v]) + smis1a[v] * (p_x[v]);
      mis[h * n_visible + v] = weighted_sum;
    }
  }
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, non-affine base
    double *restrict mish = mis + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      mish[v] /= entropy_x[v];
  }
}

static inline void tk_compressor_update_marginals (
  uint64_t cardinality,
  uint64_t *restrict bits,
  double *restrict log_p_y,
  double *restrict log_marg,
  double *restrict pseudo_counts,
  double *restrict p_y_given_x,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  tk_compressor_calculate_p_y(log_p_y, p_y_given_x, n_samples, n_hidden);
  tk_compressor_calculate_p_y_xi(cardinality, bits, log_marg, pseudo_counts, log_p_y, p_y_given_x, n_samples, n_visible, n_hidden);
}

static inline void tk_compressor_update_alpha (
  double tmin,
  double ttc,
  double lam,
  double *restrict alpha,
  double *restrict pseudo_counts,
  double *restrict tcs,
  double *restrict mis,
  double *restrict maxmis,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  for (unsigned int i = 0; i < n_hidden; i ++)
    tcs[i] = fabs(tcs[i]) * ttc + tmin;
  for (unsigned int v = 0; v < n_visible; v ++) { // not vectorized, unsupported outer form
    double max_val = 0.0;
    for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, costings not worth while
      double candidate = mis[h * n_visible + v];
      if (candidate > max_val)
        max_val = candidate;
    }
    maxmis[v] = max_val;
  }
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, non-affine base or splitting at boundary
    double t = tcs[h];
    double *restrict mish = mis + h * n_visible;
    double *restrict ps = pseudo_counts + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      ps[v] = exp(t * (mish[v] - maxmis[v]));
  }
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, non-affine
    double *restrict alphah = alpha + h * n_visible;
    double *restrict ps = pseudo_counts + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      alphah[v] = (1.0 - lam) * alphah[v] + lam * ps[v];
  }
}

static inline void tk_compressor_update_tc (
  unsigned int i,
  double *restrict log_z, // mis
  double *restrict tcs,
  double *restrict tc_history,
  unsigned int n_samples,
  unsigned int n_hidden
) {
  double sum0 = 0.0;
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, unsupported outer form
    const double *restrict lz = log_z + h * n_samples;
    double sum1 = 0.0;
    for (unsigned int s = 0; s < n_samples; s ++)
      sum1 += lz[s];
    double tc = sum1 / (double) n_samples;
    tcs[h] = tc;
    sum0 += tc;
  }
  tc_history[i] = sum0;
}

static inline void tk_compressor_calculate_latent(
  uint64_t cardinality,
  uint64_t *bits,
  double *restrict alpha,
  double *restrict p_y_given_x,
  double *restrict log_z, // mis
  double *restrict log_p_y,
  double *restrict log_p_y_given_x_unnorm,
  double *restrict log_marg,
  double *restrict sums,
  double *restrict sumsa,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  double *restrict sums0 = sums + 0 * n_hidden * n_samples;
  double *restrict sums1 = sums + 1 * n_hidden * n_samples;
  double *restrict sumsa0 = sumsa + 0 * n_hidden;
  double *restrict sumsa1 = sumsa + 1 * n_hidden;
  double *restrict lm00 = log_marg + 0 * n_hidden * n_visible;
  double *restrict lm01 = log_marg + 1 * n_hidden * n_visible;
  double *restrict lm10 = log_marg + 2 * n_hidden * n_visible;
  double *restrict lm11 = log_marg + 3 * n_hidden * n_visible;
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, unsupported outer form
    double s0 = 0.0, s1 = 0.0;
    double *restrict lm00a = lm00 + h * n_visible;
    double *restrict lm01a = lm01 + h * n_visible;
    double *restrict aph = alpha + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      s0 += aph[v] * lm00a[v];
      s1 += aph[v] * lm01a[v];
    }
    sumsa0[h] = s0;
    sumsa1[h] = s1;
  }
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, not affine
    double s0 = sumsa0[h];
    double s1 = sumsa1[h];
    double *restrict sums0a = sums0 + h * n_samples;
    double *restrict sums1a = sums1 + h * n_samples;
    for (unsigned int i = 0; i < n_samples; i ++) {
      sums0a[i] = s0;
      sums1a[i] = s1;
    }
  }
  for (unsigned int b = 0; b < cardinality; b ++) {
    unsigned int c = bits[b];
    unsigned int i = c / n_visible;
    unsigned int v = c % n_visible;
    for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, gather/scatter
      double *restrict sums0h = sums0 + h * n_samples;
      double *restrict sums1h = sums1 + h * n_samples;
      double *restrict lm00a = lm00 + h * n_visible;
      double *restrict lm01a = lm01 + h * n_visible;
      double *restrict lm10a = lm10 + h * n_visible;
      double *restrict lm11a = lm11 + h * n_visible;
      double *restrict aph = alpha + h * n_visible;
      sums0h[i] = sums0h[i] - aph[v] * lm00a[v] + aph[v] * lm10a[v]; // gather/scatter
      sums1h[i] = sums1h[i] - aph[v] * lm01a[v] + aph[v] * lm11a[v]; // gather/scatter
    }
  }
  for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, not affine
    double lpy_d0 = log_p_y[h * 2 + 0];
    double lpy_d1 = log_p_y[h * 2 + 1];
    double *restrict sums0h = sums0 + h * n_samples;
    double *restrict sums1h = sums1 + h * n_samples;
    double *restrict lpyx0 = log_p_y_given_x_unnorm + 0 * n_hidden * n_samples + h * n_samples;
    double *restrict lpyx1 = log_p_y_given_x_unnorm + 1 * n_hidden * n_samples + h * n_samples;
    for (unsigned int i = 0; i < n_samples; i ++) {
      lpyx0[i] = sums0h[i] + lpy_d0;
      lpyx1[i] = sums1h[i] + lpy_d1;
    }
  }
  tk_compressor_normalize_latent(log_z, p_y_given_x, log_p_y_given_x_unnorm, n_hidden, n_samples);
}

static inline int tk_bitmap_compress (lua_State *L)
{
  tk_compressor_t *compressor = peek_compressor(L, lua_upvalueindex(1));
  roaring64_bitmap_t *bm = peek(L, 1);
  uint64_t cardinality = roaring64_bitmap_get_cardinality(bm);
  uint64_t *bits = malloc(cardinality * sizeof(uint64_t));
  roaring64_bitmap_to_uint64_array(bm, bits);
  unsigned int n_samples = tk_lua_optunsigned(L, 2, 1);
  compressor->mis = realloc(compressor->mis, compressor->n_hidden * n_samples * sizeof(double));
  unsigned int len_sums =
    (2 * compressor->n_hidden * (n_samples < compressor->n_visible
      ? compressor->n_visible
      : n_samples));
  compressor->sums = realloc(compressor->sums, len_sums * sizeof(double));
  tk_compressor_calculate_latent(
    cardinality, bits,
    compressor->alpha,
    compressor->p_y_given_x,
    compressor->mis,
    compressor->log_p_y,
    compressor->log_p_y_given_x_unnorm,
    compressor->log_marg,
    compressor->sums,
    compressor->sumsa,
    n_samples,
    compressor->n_visible,
    compressor->n_hidden);
  free(bits);
  roaring64_bitmap_t *bm0 = roaring64_bitmap_create();
  if (bm0 == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bm0p = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *));
  *bm0p = bm0;
  luaL_getmetatable(L, MT);
  lua_setmetatable(L, -2);
  for (unsigned int h = 0; h < compressor->n_hidden; h ++) {
    for (unsigned int s = 0; s < n_samples; s ++) {
      double py0 = compressor->p_y_given_x[0 * compressor->n_hidden * n_samples + h * n_samples + s];
      double py1 = compressor->p_y_given_x[1 * compressor->n_hidden * n_samples + h * n_samples + s];
      if (py1 > py0)
        roaring64_bitmap_add(bm0, s * compressor->n_hidden + h);
    }
  }
  return 1;
}

static uint64_t const multiplier = 6364136223846793005u;
__thread uint64_t mcg_state = 0xcafef00dd15ea5e5u;

static inline uint32_t fast_rand ()
{
  uint64_t x = mcg_state;
  unsigned int count = (unsigned int) (x >> 61);
  mcg_state = x * multiplier;
  return (uint32_t) ((x ^ x >> 22) >> (22 + count));
}

static inline double fast_drand ()
{
  return ((double)fast_rand()) / ((double)UINT32_MAX);
}

static inline bool tk_compressor_converged (
  unsigned int i,
  double eps,
  double *tc_history
) {
  if (i < 10)
    return false;
  double m0 = 0.0;
  for (unsigned int j = i - 10; j < i - 5; j ++)
    m0 += tc_history[j];
  m0 = m0 / 5;
  double m1 = 0.0;
  for (unsigned int j = i - 5; j < i; j ++)
    m1 += tc_history[j];
  m1 = m1 / 5;
  return fabs(-m0 + m1) < eps;
}

static inline void tk_compressor_init_alpha (
  double *alpha,
  unsigned int n_hidden,
  unsigned int n_visible
) {
  for (unsigned int i = 0; i < n_hidden * n_visible; ++i)
    alpha[i] = 0.5 + 0.5 * fast_drand();
}

static inline void tk_compressor_init_tcs (
  double *tcs,
  unsigned int n_hidden
) {
  for (unsigned int i = 0; i < n_hidden; i ++)
    tcs[i] = 0.0;
}

static inline void tk_compressor_init_log_p_y_given_x_unnorm (
  double *log_p_y_given_x_unnorm,
  unsigned int n_hidden,
  unsigned int n_samples
) {
  double log_dim_hidden = -log(2);
  for (unsigned int i = 0; i < 2 * n_hidden * n_samples; ++i)
    log_p_y_given_x_unnorm[i] = log_dim_hidden * (0.5 + fast_drand());
}

static inline void tk_compressor_init (
  lua_State *L,
  tk_compressor_t *compressor,
  roaring64_bitmap_t *bm,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int max_iter,
  double eps,
  int i_each
) {
  compressor->n_visible = n_visible;
  compressor->n_hidden = n_hidden;
  compressor->lam = 0.3;
  compressor->tmin = 1.0;
  compressor->ttc = 500.0;
  compressor->tcs = malloc(compressor->n_hidden * sizeof(double));
  compressor->tc_history = malloc(max_iter * sizeof(double));
  compressor->alpha = malloc(compressor->n_hidden * compressor->n_visible * sizeof(double));
  compressor->p_y_given_x = malloc(2 * compressor->n_hidden * n_samples * sizeof(double));
  compressor->log_p_y = malloc(2 * compressor->n_hidden * sizeof(double));
  compressor->log_p_y_given_x_unnorm = malloc(2 * compressor->n_hidden * n_samples * sizeof(double));
  compressor->p_x = malloc(compressor->n_visible * sizeof(double));
  compressor->entropy_x = malloc(compressor->n_visible * sizeof(double));
  compressor->log_marg = malloc(2 * 2 * compressor->n_hidden * compressor->n_visible * sizeof(double));
  compressor->pseudo_counts = malloc(2 * 2 * compressor->n_hidden * compressor->n_visible * sizeof(double));
  unsigned int len_mis = compressor->n_hidden * compressor->n_visible;
  if (len_mis < (compressor->n_hidden * n_samples))
    len_mis = compressor->n_hidden * n_samples;
  compressor->mis = malloc(len_mis * sizeof(double));
  compressor->smis = malloc(compressor->n_hidden * compressor->n_visible * 2 * sizeof(double));
  compressor->maxmis = malloc(compressor->n_visible * sizeof(double));
  compressor->sums = malloc(2 * n_samples * n_hidden * sizeof(double));
  compressor->sumsa = malloc(2 * n_hidden * sizeof(double));
  uint64_t cardinality = roaring64_bitmap_get_cardinality(bm);
  uint64_t *bits = malloc(cardinality * sizeof(uint64_t));
  roaring64_bitmap_to_uint64_array(bm, bits);
  tk_compressor_data_stats(
    cardinality, bits,
    compressor->p_x,
    compressor->entropy_x,
    n_samples,
    compressor->n_visible);
  tk_compressor_init_tcs(
    compressor->tcs,
    compressor->n_hidden);
  tk_compressor_init_alpha(
    compressor->alpha,
    compressor->n_hidden,
    compressor->n_visible);
  tk_compressor_init_log_p_y_given_x_unnorm(
    compressor->log_p_y_given_x_unnorm,
    compressor->n_hidden,
    n_samples);
  tk_compressor_normalize_latent(
    compressor->mis,
    compressor->p_y_given_x,
    compressor->log_p_y_given_x_unnorm,
    compressor->n_hidden,
    n_samples);
  unsigned int i = 0;
  while (i < max_iter) {
    tk_compressor_update_marginals(
      cardinality, bits,
      compressor->log_p_y,
      compressor->log_marg,
      compressor->pseudo_counts,
      compressor->p_y_given_x,
      n_samples,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_calculate_mis(
      compressor->mis,
      compressor->smis,
      compressor->pseudo_counts,
      compressor->log_p_y,
      compressor->log_marg,
      compressor->p_x,
      compressor->entropy_x,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_update_alpha(
      compressor->tmin,
      compressor->ttc,
      compressor->lam,
      compressor->alpha,
      compressor->pseudo_counts,
      compressor->tcs,
      compressor->mis,
      compressor->maxmis,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_calculate_latent(
      cardinality, bits,
      compressor->alpha,
      compressor->p_y_given_x,
      compressor->mis,
      compressor->log_p_y,
      compressor->log_p_y_given_x_unnorm,
      compressor->log_marg,
      compressor->sums,
      compressor->sumsa,
      n_samples,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_update_tc(
      i,
      compressor->mis,
      compressor->tcs,
      compressor->tc_history,
      n_samples,
      compressor->n_hidden);
    bool converged = tk_compressor_converged(i, eps, compressor->tc_history);
    if (i_each > -1) {
      lua_pushvalue(L, i_each);
      lua_pushinteger(L, i + 1);
      lua_pushnumber(L, compressor->tc_history[i]);
      lua_pushboolean(L, converged);
      lua_call(L, 3, 0);
    }
    if (converged)
      break;
    else
      i ++;
  }
  tk_compressor_shrink(compressor);
  free(bits);
}

static inline int tk_bitmap_compressor (lua_State *L)
{
  roaring64_bitmap_t *bm = peek(L, 1);
  unsigned int n_samples = tk_lua_checkunsigned(L, 2);
  unsigned int n_visible = tk_lua_checkunsigned(L, 3);
  unsigned int n_hidden = tk_lua_checkunsigned(L, 4);
  unsigned int max_iter = tk_lua_checkunsigned(L, 5);
  double eps = luaL_optnumber(L, 6, 1e-6);
  unsigned int i_each = -1;
  if (lua_type(L, 7) != LUA_TNIL)
    i_each = tk_lua_absindex(L, i_each);
  if (max_iter < 10)
    max_iter = 10;
  tk_compressor_t *compressor = (tk_compressor_t *)
    lua_newuserdata(L, sizeof(tk_compressor_t));
  luaL_getmetatable(L, MT_COMPRESSOR); // t, n, b, mt
  lua_setmetatable(L, -2); // t, n, b
  tk_compressor_init(L, compressor, bm, n_samples, n_visible, n_hidden, max_iter, eps, i_each);
  lua_pushcclosure(L, tk_bitmap_compress, 1);
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
  { "extend", tk_bitmap_extend },
  { "bits", tk_bitmap_bits },
  { "compressor", tk_bitmap_compressor },
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
  luaL_newmetatable(L, MT_COMPRESSOR); // t mt
  lua_pushcfunction(L, tk_compressor_destroy); // t mt fn
  lua_setfield(L, -2, "__gc"); // t mt
  lua_pop(L, 1); // t
  return 1;
}
