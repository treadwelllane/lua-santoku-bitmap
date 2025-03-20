#include "lua.h"
#include "lauxlib.h"
#include "roaring.h"
#include "roaring.c"
#include <string.h>
#include <math.h>
#include <float.h>

#define MT "santoku_bitmap"

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

static void tk_bitmap_iterate (
  roaring64_bitmap_t *bm,
  bool (*it)(uint64_t, void *),
  void *udata
) {
  for (uint64_t i = roaring64_bitmap_minimum(bm); i <= roaring64_bitmap_maximum(bm); i ++)
    if (roaring64_bitmap_contains(bm, i))
      if (!it(i, udata))
        break;
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
    bits = luaL_checkinteger(L, 2);
    if (bits < 0) {
      return luaL_error(L, "number of bits can't be negative");
    } else if (bits > UINT64_MAX) {
      return luaL_error(L, "number of bits can't be greater than UINT64_MAX");
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
  tk_bitmap_iterate(bm, tk_bitmap_raw_iter, &state);
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
  tk_bitmap_iterate(bm, tk_bitmap_bits_iter, L); // bm t bits n
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
  tk_bitmap_iterate(bm1, tk_bitmap_extend_iter, &state);
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
  double *alphaopt;
  double *log_marg;
  double *log_marg_xi;
  double *log_p_y;
  double *log_p_y_given_x_unnorm;
  double *log_z;
  double *maxmis;
  double *mis;
  double *p_x;
  double *entropy_x;
  double *p_y_given_x;
  double *pseudo_counts;
  double *smis;
  double *tc_history;
  double *tcs;
  double *vec;
  double log_dim_hidden;
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
  if (compressor->alphaopt) {
    free(compressor->alphaopt);
    compressor->alphaopt = NULL;
  }
  if (compressor->entropy_x) {
    free(compressor->entropy_x);
    compressor->entropy_x = NULL;
  }
  if (compressor->log_marg_xi) {
    free(compressor->log_marg_xi);
    compressor->log_marg_xi = NULL;
  }
  if (compressor->maxmis) {
    free(compressor->maxmis);
    compressor->maxmis = NULL;
  }
  if (compressor->mis) {
    free(compressor->mis);
    compressor->mis = NULL;
  }
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
  if (compressor->vec) {
    free(compressor->vec);
    compressor->vec = NULL;
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
  if (compressor->log_z) {
    free(compressor->log_z);
    compressor->log_z = NULL;
  }
  if (compressor->p_y_given_x) {
    free(compressor->p_y_given_x);
    compressor->p_y_given_x = NULL;
  }
  return 1;
}

static inline int tk_bitmap_compress (lua_State *);

static inline void tk_compressor_logsumexp (
  double *log_p_y_given_x_unnorm,
  unsigned int n_hidden,
  unsigned int n_samples,
  double *log_z
) {
  for (unsigned int s = 0; s < n_samples; s++) {
    for (unsigned int h = 0; h < n_hidden; h++) {
      double a = log_p_y_given_x_unnorm[h * n_samples * 2 + s * 2 + 0];
      double b = log_p_y_given_x_unnorm[h * n_samples * 2 + s * 2 + 1];
      double max_ab  = (a > b)? a : b;
      double sum_exp = exp(a - max_ab) + exp(b - max_ab);
      log_z[h * n_samples + s] = max_ab + log(sum_exp);
    }
  }
}

static void tk_compressor_normalize_latent (
  double *log_z,
  double *p_y_given_x,
  double *log_p_y_given_x_unnorm,
  unsigned int n_hidden,
  unsigned int n_samples
) {
  tk_compressor_logsumexp(log_p_y_given_x_unnorm, n_hidden, n_samples, log_z);
  for (int i = 0; i < n_hidden; ++i) {
    for (int j = 0; j < n_samples; ++j) {
      for (int k = 0; k < 2; ++k) {
        p_y_given_x[i * n_samples  * 2 + j * 2 + k] =
          exp(log_p_y_given_x_unnorm[i * n_samples * 2 + j * 2 + k] -
              log_z[i * n_samples + j]);
      }
    }
  }
}

static inline void tk_compressor_data_stats (
  roaring64_bitmap_t *bm,
  double *p_x,
  double *entropy_x,
  double eps,
  unsigned int n_samples,
  unsigned int n_visible
) {
  for (unsigned int i = 0; i < n_visible; i++) {
    double count1 = 0.0;
    for (unsigned int j = 0; j < n_samples; j++)
      if (roaring64_bitmap_contains(bm, j * n_visible + i))
        count1 += 1.0;
    double count0 = (double) n_samples - count1;
    double total = count0 + count1;
    p_x[i * 2 + 0] = count0 / total;
    p_x[i * 2 + 1] = count1 / total;
  }
  for (unsigned int i = 0; i < n_visible; i++) {
    double p0 = p_x[i * 2 + 0];
    double p1 = p_x[i * 2 + 1];
    double entropy = 0.0;
    if (p0 > 0.0)
      entropy -= p0 * log(p0);
    if (p1 > 0.0)
      entropy -= p1 * log(p1);
    entropy_x[i] = (entropy > 0.0) ? entropy : 1e-10;
  }
}

static inline void tk_compressor_calculate_p_y (
  double *log_p_y,
  double *p_y_given_x,
  unsigned int n_samples,
  unsigned int n_hidden
) {
  for (unsigned int h = 0; h < n_hidden; ++h) {
    double pseudo_counts_0 = 0.001;
    double pseudo_counts_1 = 0.001;
    for (unsigned int s = 0; s < n_samples; ++s) {
      pseudo_counts_0 += p_y_given_x[(h * n_samples + s) * 2 + 0];
      pseudo_counts_1 += p_y_given_x[(h * n_samples + s) * 2 + 1];
    }
    double sum_pseudo_counts = pseudo_counts_0 + pseudo_counts_1;
    log_p_y[h * 2 + 0] = log(pseudo_counts_0) - log(sum_pseudo_counts);
    log_p_y[h * 2 + 1] = log(pseudo_counts_1) - log(sum_pseudo_counts);
  }
}

static void tk_compressor_calculate_p_y_xi (
  roaring64_bitmap_t *bm,
  double *log_marg,
  double *pseudo_counts,
  double *p_y_given_x,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  const unsigned int dim_hidden = 2;
  const unsigned int n_events = n_visible * 2;
  const double smoothing = 0.001;
  for (unsigned int j = 0; j < n_visible; j++) {
    unsigned int event0 = 2 * j;
    unsigned int event1 = event0 + 1;
    for (unsigned int h = 0; h < n_hidden; h++) {
      double sum0_dim0 = 0.0, sum0_dim1 = 0.0;
      double sum1_dim0 = 0.0, sum1_dim1 = 0.0;
      for (unsigned int k = 0; k < n_samples; k++) {
        unsigned int base_idx = (h * n_samples + k) * dim_hidden;
        double py0 = p_y_given_x[base_idx + 0];
        double py1 = p_y_given_x[base_idx + 1];
        uint64_t index = ((uint64_t)k) * n_visible + j;
        if (roaring64_bitmap_contains(bm, index)) {
          sum1_dim0 += py0;
          sum1_dim1 += py1;
        } else {
          sum0_dim0 += py0;
          sum0_dim1 += py1;
        }
      }
      unsigned int idx0 = (h * n_events + event0) * dim_hidden;
      unsigned int idx1 = (h * n_events + event1) * dim_hidden;
      pseudo_counts[idx0 + 0] = sum0_dim0 + smoothing;
      pseudo_counts[idx0 + 1] = sum0_dim1 + smoothing;
      pseudo_counts[idx1 + 0] = sum1_dim0 + smoothing;
      pseudo_counts[idx1 + 1] = sum1_dim1 + smoothing;
    }
  }
  for (unsigned int h = 0; h < n_hidden; h++) {
    for (unsigned int e = 0; e < n_events; e++) {
      unsigned int base = (h * n_events + e) * dim_hidden;
      double total = pseudo_counts[base] + pseudo_counts[base + 1];
      double log_total = log(total);
      for (unsigned int d = 0; d < dim_hidden; d++) {
        unsigned int idx = base + d;
        log_marg[idx] = log(pseudo_counts[idx]) - log_total;
      }
    }
  }
}

static inline void tk_compressor_calculate_mis (
  double *mis,
  double *smis,
  double *vec,
  double *log_p_y,
  double *log_marg,
  double *p_x,
  double *entropy_x,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  const unsigned int n_events = n_visible * 2;
  for (unsigned int h = 0; h < n_hidden; h++) {
    for (unsigned int e = 0; e < n_events; e++) {
      for (unsigned int d = 0; d < 2; d++) {
        unsigned int idx = (h * n_events + e) * 2 + d;
        vec[idx] = exp(log_marg[idx] + log_p_y[h * 2 + d]);
      }
    }
  }
  for (unsigned int h = 0; h < n_hidden; h++) {
    for (unsigned int e = 0; e < n_events; e++) {
      double sum = 0.0;
      for (unsigned int d = 0; d < 2; d++) {
        unsigned int idx = (h * n_events + e) * 2 + d;
        sum += vec[idx] * log_marg[idx];
      }
      smis[h * n_events + e] = sum;
    }
  }
  for (unsigned int h = 0; h < n_hidden; h++) {
    for (unsigned int j = 0; j < n_visible; j++) {
      unsigned int idx0 = h * n_events + j * 2 + 0;
      unsigned int idx1 = h * n_events + j * 2 + 1;
      double weighted_sum =
        smis[idx0] * p_x[j * 2 + 0] +
        smis[idx1] * p_x[j * 2 + 1];
      mis[h * n_visible + j] = weighted_sum;
    }
  }
  for (unsigned int h = 0; h < n_hidden; h++) {
    for (unsigned int j = 0; j < n_visible; j++) {
      mis[h * n_visible + j] /= entropy_x[j];
    }
  }
}

static inline void tk_compressor_update_marginals (
  roaring64_bitmap_t *bm,
  double *log_p_y,
  double *log_marg,
  double *pseudo_counts,
  double *log_marg_xi,
  double *p_y_given_x,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  tk_compressor_calculate_p_y(log_p_y, p_y_given_x, n_samples, n_hidden);
  tk_compressor_calculate_p_y_xi(bm, log_marg, pseudo_counts, p_y_given_x, n_samples, n_visible, n_hidden);
  for (unsigned int h = 0; h < n_hidden; h++) {
    for (unsigned int e = 0; e < n_visible * 2; e++) {
      for (unsigned int d = 0; d < 2; d++) {
        unsigned int idx = (h * n_visible * 2 + e) * 2 + d;
        log_marg[idx] -= log_p_y[h * 2 + d];
      }
    }
  }
}

static inline void tk_compressor_update_alpha (
  double tmin,
  double ttc,
  double lam,
  double *alpha,
  double *alphaopt,
  double *tcs,
  double *mis,
  double *maxmis,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  for (unsigned int i = 0; i < n_hidden; i ++)
    tcs[i] = fabs(tcs[i]) * ttc + tmin;
  for (unsigned int v = 0; v < n_visible; v ++) {
    maxmis[v] = 0;
    for (unsigned int h = 0; h < n_hidden; h ++)
      if (mis[h * n_visible + v] > maxmis[v])
        maxmis[v] = mis[h * n_visible + v];
  }
  for (int i = 0; i < n_hidden; i ++)
    for (int j = 0; j < n_visible; j ++)
      alphaopt[i * n_visible + j] = exp(tcs[i] * (mis[i * n_visible + j] - maxmis[j]));
  for (int i = 0; i < n_hidden; i ++)
    for (int j = 0; j < n_visible; j ++)
      alpha[i * n_visible + j] = (1.0 - lam) * alpha[i * n_visible + j] + lam * alphaopt[i * n_visible + j];
}

static inline void tk_compressor_update_tc (
  unsigned int i,
  double *log_z,
  double *tcs,
  double *tc_history,
  unsigned int n_samples,
  unsigned int n_hidden
) {
  double sum0 = 0.0;
  for (unsigned int h = 0; h < n_hidden; h ++) {
    double sum1 = 0.0;
    for (unsigned int s = 0; s < n_samples; s ++)
      sum1 += log_z[h * n_samples + s];
    tcs[h] = sum1 / (double) n_samples;
    sum0 += tcs[h];
  }
  tc_history[i] = sum0;
}

static inline void tk_compressor_calculate_latent (
  roaring64_bitmap_t *bm,
  double *alpha,
  double *p_y_given_x,
  double *log_z,
  double *log_p_y,
  double *log_p_y_given_x_unnorm,
  double *log_marg,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden
) {
  unsigned int n_events = n_visible * 2;
  for (unsigned int i = 0; i < n_samples; i++) {
    for (unsigned int h = 0; h < n_hidden; h++) {
      for (unsigned int d = 0; d < 2; d++) {
        double s = 0.0;
        for (unsigned int e = 0; e < n_events; e++) {
          unsigned int j = e / 2;
          unsigned int state = e % 2;
          int has_feature = roaring64_bitmap_contains(bm, ((uint64_t)i) * n_visible + j);
          double X_val = state ? (has_feature ? 1 : 0) : (has_feature ? 0 : 1);
          s += X_val * alpha[h * n_visible + j] * log_marg[h * n_events * 2 + e * 2 + d];
        }
        log_p_y_given_x_unnorm[h * n_samples * 2 + i * 2 + d] = s + log_p_y[h * 2 + d];
      }
    }
  }
  tk_compressor_normalize_latent(log_z, p_y_given_x, log_p_y_given_x_unnorm, n_hidden, n_samples);
}

static inline int tk_bitmap_compress (lua_State *L)
{
  tk_compressor_t *compressor = peek_compressor(L, lua_upvalueindex(1));
  roaring64_bitmap_t *bm = peek(L, 1);
  unsigned int n_samples = tk_lua_optunsigned(L, 2, 1);
  tk_compressor_calculate_latent(
    bm,
    compressor->alpha,
    compressor->p_y_given_x,
    compressor->log_z,
    compressor->log_p_y,
    compressor->log_p_y_given_x_unnorm,
    compressor->log_marg,
    n_samples,
    compressor->n_visible,
    compressor->n_hidden);
  roaring64_bitmap_t *bm0 = roaring64_bitmap_create();
  if (bm0 == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bm0p = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *));
  *bm0p = bm0;
  luaL_getmetatable(L, MT);
  lua_setmetatable(L, -2);
  for (unsigned int h = 0; h < compressor->n_hidden; ++h) {
    for (unsigned int s = 0; s < n_samples; ++s) {
      double py0 = compressor->p_y_given_x[h * n_samples * 2 + s * 2 + 0];
      double py1 = compressor->p_y_given_x[h * n_samples * 2 + s * 2 + 1];
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

static inline void tk_compressor_init (
  lua_State *L,
  tk_compressor_t *compressor,
  roaring64_bitmap_t *bm,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int max_iter
) {
  compressor->n_visible = n_visible;
  compressor->n_hidden = n_hidden;
  compressor->eps = 1e-6;
  compressor->lam = 0.3;
  compressor->tmin = 1.0;
  compressor->ttc = 500.0;
  compressor->tcs = malloc(compressor->n_hidden * sizeof(double));
  for (unsigned int i = 0; i < compressor->n_hidden; i ++)
    compressor->tcs[i] = 0.0;
  compressor->tc_history = malloc(max_iter * sizeof(double));
  compressor->alpha = malloc(compressor->n_hidden * compressor->n_visible * sizeof(double));
  compressor->alphaopt = malloc(compressor->n_hidden * compressor->n_visible * sizeof(double));
  for (int i = 0; i < compressor->n_hidden; ++i)
    for (int j = 0; j < compressor->n_visible; ++j)
      compressor->alpha[i * compressor->n_visible + j] = 0.5 + 0.5 * fast_drand();
  double log_dim_hidden = -log(2);
  compressor->log_p_y_given_x_unnorm = malloc(compressor->n_hidden * n_samples * 2 * sizeof(double));
  for (int i = 0; i < compressor->n_hidden; ++i)
    for (int j = 0; j < n_samples; ++j)
      for (int k = 0; k < 2; ++k)
        compressor->log_p_y_given_x_unnorm[i * n_samples * 2 + j * 2 + k] = log_dim_hidden * (0.5 + fast_drand());
  compressor->log_z = malloc(compressor->n_hidden * n_samples * sizeof(double));
  compressor->log_p_y = malloc(compressor->n_hidden * 2 * sizeof(double));
  compressor->p_y_given_x = malloc(compressor->n_hidden * n_samples * 2 * sizeof(double));
  tk_compressor_normalize_latent(
    compressor->log_z,
    compressor->p_y_given_x,
    compressor->log_p_y_given_x_unnorm,
    compressor->n_hidden,
    n_samples);
  compressor->p_x = malloc(compressor->n_visible * 2 * sizeof(double));
  compressor->entropy_x = malloc(compressor->n_visible * sizeof(double));
  tk_compressor_data_stats(
    bm,
    compressor->p_x,
    compressor->entropy_x,
    compressor->eps,
    n_samples,
    compressor->n_visible);
  compressor->log_marg = malloc(compressor->n_hidden * compressor->n_visible * 2 * 2 * sizeof(double));
  compressor->vec = malloc(compressor->n_hidden * compressor->n_visible * 2 * 2 * sizeof(double));
  compressor->log_marg_xi = malloc(compressor->n_hidden * compressor->n_visible * 2 * sizeof(double));
  compressor->pseudo_counts = malloc(compressor->n_hidden * compressor->n_visible * 2 * 2 * sizeof(double));
  compressor->mis = malloc(compressor->n_hidden * compressor->n_visible * sizeof(double));
  compressor->smis = malloc(compressor->n_hidden * compressor->n_visible * 2 * sizeof(double));
  compressor->maxmis = malloc(compressor->n_visible * sizeof(double));
  unsigned int i = 0;
  while (i < max_iter) {
    tk_compressor_update_marginals(
      bm,
      compressor->log_p_y,
      compressor->log_marg,
      compressor->pseudo_counts,
      compressor->log_marg_xi,
      compressor->p_y_given_x,
      n_samples,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_calculate_mis(
      compressor->mis,
      compressor->smis,
      compressor->vec,
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
      compressor->alphaopt,
      compressor->tcs,
      compressor->mis,
      compressor->maxmis,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_calculate_latent(
      bm,
      compressor->alpha,
      compressor->p_y_given_x,
      compressor->log_z,
      compressor->log_p_y,
      compressor->log_p_y_given_x_unnorm,
      compressor->log_marg,
      n_samples,
      compressor->n_visible,
      compressor->n_hidden);
    tk_compressor_update_tc(
      i,
      compressor->log_z,
      compressor->tcs,
      compressor->tc_history,
      n_samples,
      compressor->n_hidden);
    if (tk_compressor_converged(i, compressor->eps, compressor->tc_history))
      break;
    else
      i ++;
  }
  tk_compressor_shrink(compressor);
}

static inline int tk_bitmap_compressor (lua_State *L)
{
  lua_pushnil(L);
  roaring64_bitmap_t *bm = peek(L, 1);
  unsigned int n_rows = tk_lua_checkunsigned(L, 2);
  unsigned int n_cols = tk_lua_checkunsigned(L, 3);
  unsigned int n_reduced = tk_lua_checkunsigned(L, 4);
  unsigned int max_iter = tk_lua_checkunsigned(L, 5);
  if (max_iter < 10)
    max_iter = 10;
  tk_compressor_t *compressor = (tk_compressor_t *)
    lua_newuserdata(L, sizeof(tk_compressor_t));
  luaL_getmetatable(L, MT_COMPRESSOR); // t, n, b, mt
  lua_setmetatable(L, -2); // t, n, b
  tk_compressor_init(L, compressor, bm, n_rows, n_cols, n_reduced, max_iter);
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
