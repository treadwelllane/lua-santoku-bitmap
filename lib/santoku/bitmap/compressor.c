#include "lua.h"
#include "lauxlib.h"
#include "khash.h"
#include "roaring.h"
#include "roaring.c"
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define MT_BITMAP "santoku_bitmap"
#define MT_COMPRESSOR "santoku_bitmap_compressor"

typedef struct {
  unsigned int d0;
  unsigned int d1;
  unsigned int h;
  unsigned int v;
} mtx_key_t;

static inline int mtx_hash (mtx_key_t k)
{
  khint_t a = kh_int_hash_func(k.d0);
  khint_t b = kh_int_hash_func(k.d1);
  khint_t c = kh_int_hash_func(k.h);
  khint_t d = kh_int_hash_func(k.v);
  khint_t hash = a;
  hash = a * 31 + b;
  hash = b * 31 + c;
  hash = c * 31 + d;
  return hash;
}

static inline int mtx_equals (mtx_key_t a, mtx_key_t b)
{
  return a.d0 == b.d0 && a.d1 == b.d1 && a.h == b.h && a.v == b.v;
}

KHASH_INIT(mtx, mtx_key_t, double, 1, mtx_hash, mtx_equals);
typedef khash_t(mtx) tb_mtx_t;

static inline void mtx_add (tb_mtx_t *m, unsigned int d0, unsigned int d1, unsigned int h, unsigned int v, double x)
{
  mtx_key_t i = { .d0 = d0, .d1 = d1, .h = h, .v = v };
  int absent;
  khint_t k = kh_put(mtx, m, i, &absent);
  kh_value(m, k) = absent ? x : kh_value(m, k) + x;
  if (fabs(kh_value(m, k)) < 1e-8)
    kh_del(mtx, m, k);
}

static inline void mtx_set (tb_mtx_t *m, unsigned int d0, unsigned int d1, unsigned int h, unsigned int v, double x)
{
  mtx_key_t i = { .d0 = d0, .d1 = d1, .h = h, .v = v };
  if (fabs(x) < 1e-8) {
    khint_t k = kh_get(mtx, m, i);
    if (k != kh_end(m))
      kh_del(mtx, m, k);
  } else {
    int absent;
    khint_t k = kh_put(mtx, m, i, &absent);
    kh_value(m, k) = x;
  }
}

static inline double mtx_get (tb_mtx_t *m, unsigned int d0, unsigned int d1, unsigned int h, unsigned int v)
{
  mtx_key_t i = { .d0 = d0, .d1 = d1, .h = h, .v = v };
  khint_t k = kh_get(mtx, m, i);
  return k == kh_end(m) ? 0 : kh_value(m, k);
}

#define lm_clear(C) kh_clear(mtx, C->lm)
#define lm_add(C, d0, d1, h, v, x) mtx_add((C)->lm, d0, d1, h, v, x)
#define lm_set(C, d0, d1, h, v, x) mtx_set((C)->lm, d0, d1, h, v, x)
#define lm_del(C, d0, d1, h, v) mtx_set((C)->lm, d0, d1, h, v, 0)
#define lm_get(C, d0, d1, h, v) mtx_get((C)->lm, d0, d1, h, v)

#define pc_clear(C) kh_clear(mtx, C->cnt)
#define pc_add(C, d0, d1, h, v, x) mtx_add((C)->cnt, d0, d1, h, v, x)
#define pc_set(C, d0, d1, h, v, x) mtx_set((C)->cnt, d0, d1, h, v, x)
#define pc_del(C, d0, d1, h, v) mtx_set((C)->cnt, d0, d1, h, v, 0)
#define pc_get(C, d0, d1, h, v) mtx_get((C)->cnt, d0, d1, h, v)

typedef enum {
  TK_CMP_INIT,
  TK_CMP_INIT_ALPHA,
  TK_CMP_INIT_TCS,
  TK_CMP_INIT_PYX_UNNORM,
  TK_CMP_MARGINALS,
  TK_CMP_MAXMIS,
  TK_CMP_ALPHA,
  TK_CMP_LATENT_SUMSA,
  TK_CMP_LATENT_SUMS,
  TK_CMP_LATENT_PY,
  TK_CMP_LATENT_NORM,
  TK_CMP_UPDATE_TC,
  TK_CMP_DONE
} tk_compressor_stage_t;

struct tk_compressor_s;
typedef struct tk_compressor_s tk_compressor_t;

typedef struct {
  tk_compressor_t *C;
  unsigned int cardinality;
  unsigned int n_samples;
  tk_compressor_stage_t stage;
  unsigned int hfirst;
  unsigned int hlast;
  unsigned int vfirst;
  unsigned int vlast;
  unsigned int bfirst;
  unsigned int blast;
} tk_compressor_thread_data_t;

typedef struct tk_compressor_s {
  bool trained; // already trained
  bool destroyed;
  double *alpha;
  tb_mtx_t *lm;
  tb_mtx_t *cnt;
  double *tmp;
  double *log_py;
  double *log_pyx_unnorm;
  double *maxmis;
  double *mis;
  double *sums;
  double *sumsa;
  uint64_t *samples;
  unsigned int *visibles;
  double *px;
  double *entropy_x;
  double *pyx;
  double *smis;
  double last_tc;
  double *tcs;
  double lam;
  double tmin;
  double ttc;
  unsigned int n_visible;
  unsigned int n_hidden;
  pthread_mutex_t mutex;
  pthread_cond_t cond_stage;
  pthread_cond_t cond_done;
  unsigned int n_threads;
  unsigned int n_threads_done;
  tk_compressor_stage_t stage;
  pthread_t *threads;
  tk_compressor_thread_data_t *thread_data;
  bool created_threads;
} tk_compressor_t;

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

static inline int tk_lua_error (lua_State *L, const char *err)
{
  lua_pushstring(L, err);
  tk_lua_callmod(L, 1, 0, "santoku.error", "error");
  return 0;
}

static inline FILE *tk_lua_fmemopen (lua_State *L, char *data, size_t size, const char *flag)
{
  FILE *fh = fmemopen(data, size, flag);
  if (fh) return fh;
  int e = errno;
  lua_settop(L, 0);
  lua_pushstring(L, "Error opening string as file");
  lua_pushstring(L, strerror(e));
  lua_pushinteger(L, e);
  tk_lua_callmod(L, 3, 0, "santoku.error", "error");
  return NULL;
}

static inline FILE *tk_lua_fopen (lua_State *L, const char *fp, const char *flag)
{
  FILE *fh = fopen(fp, flag);
  if (fh) return fh;
  int e = errno;
  lua_settop(L, 0);
  lua_pushstring(L, "Error opening file");
  lua_pushstring(L, fp);
  lua_pushstring(L, strerror(e));
  lua_pushinteger(L, e);
  tk_lua_callmod(L, 4, 0, "santoku.error", "error");
  return NULL;
}

static inline void tk_lua_fclose (lua_State *L, FILE *fh)
{
  if (!fclose(fh)) return;
  int e = errno;
  lua_settop(L, 0);
  lua_pushstring(L, "Error closing file");
  lua_pushstring(L, strerror(e));
  lua_pushinteger(L, e);
  tk_lua_callmod(L, 3, 0, "santoku.error", "error");
}

static inline void tk_lua_fwrite (lua_State *L, void *data, size_t size, size_t memb, FILE *fh)
{
  fwrite(data, size, memb, fh);
  if (!ferror(fh)) return;
  int e = errno;
  lua_settop(L, 0);
  lua_pushstring(L, "Error writing to file");
  lua_pushstring(L, strerror(e));
  lua_pushinteger(L, e);
  tk_lua_callmod(L, 3, 0, "santoku.error", "error");
}

static inline void tk_lua_fread (lua_State *L, void *data, size_t size, size_t memb, FILE *fh)
{
  size_t r = fread(data, size, memb, fh);
  if (!ferror(fh) || !r) return;
  int e = errno;
  lua_settop(L, 0);
  lua_pushstring(L, "Error reading from file");
  lua_pushstring(L, strerror(e));
  lua_pushinteger(L, e);
  tk_lua_callmod(L, 3, 0, "santoku.error", "error");
}

static inline int tk_lua_errno (lua_State *L, int err)
{
  lua_pushstring(L, strerror(errno));
  lua_pushinteger(L, err);
  tk_lua_callmod(L, 2, 0, "santoku.error", "error");
  return 0;
}

static inline int tk_lua_errmalloc (lua_State *L)
{
  lua_pushstring(L, "Error in malloc");
  tk_lua_callmod(L, 1, 0, "santoku.error", "error");
  return 0;
}

static inline FILE *tk_lua_tmpfile (lua_State *L)
{
  FILE *fh = tmpfile();
  if (fh) return fh;
  int e = errno;
  lua_settop(L, 0);
  lua_pushstring(L, "Error opening tmpfile");
  lua_pushstring(L, strerror(e));
  lua_pushinteger(L, e);
  tk_lua_callmod(L, 3, 0, "santoku.error", "error");
  return NULL;
}

static inline char *tk_lua_fslurp (lua_State *L, FILE *fh, size_t *len)
{
  if (fseek(fh, 0, SEEK_END) != 0) {
    tk_lua_errno(L, errno);
    return NULL;
  }
  long size = ftell(fh);
  if (size < 0) {
    tk_lua_errno(L, errno);
    return NULL;
  }
  if (fseek(fh, 0, SEEK_SET) != 0) {
    tk_lua_errno(L, errno);
    return NULL;
  }
  char *buffer = malloc((size_t) size);
  if (!buffer) {
    tk_lua_errmalloc(L);
    return NULL;
  }
  if (fread(buffer, 1, (size_t) size, fh) != (size_t) size) {
    free(buffer);
    tk_lua_errno(L, errno);
    return NULL;
  }
  *len = (size_t) size;
  return buffer;
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

static inline lua_Number tk_lua_foptnumber (lua_State *L, int i, char *field, double d)
{
  lua_getfield(L, i, field);
  lua_Number n = luaL_optnumber(L, -1, d);
  lua_pop(L, 1);
  return n;
}

static inline lua_Integer tk_lua_fcheckunsigned (lua_State *L, int i, char *field)
{
  lua_getfield(L, i, field);
  lua_Integer n = tk_lua_checkunsigned(L, -1);
  lua_pop(L, 1);
  return n;
}

static inline unsigned int tk_lua_optunsigned (lua_State *L, int i, unsigned int def)
{
  if (lua_type(L, i) < 1)
    return def;
  return tk_lua_checkunsigned(L, i);
}

static inline lua_Integer tk_lua_ftype (lua_State *L, int i, char *field)
{
  lua_getfield(L, i, field);
  int t = lua_type(L, -1);
  lua_pop(L, 1);
  return t;
}

static inline void *tk_lua_fcheckuserdata (lua_State *L, int i, char *field, char *mt)
{
  lua_getfield(L, i, field);
  void *p = luaL_checkudata(L, -1, mt);
  lua_pop(L, 1);
  return p;
}

static roaring64_bitmap_t *peek_bitmap (lua_State *L, int i)
{
  return *((roaring64_bitmap_t **) luaL_checkudata(L, i, MT_BITMAP));
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

static inline void *tk_realloc (
  lua_State *L,
  void *p,
  size_t s
) {
  p = realloc(p, s);
  if (!p) {
    tk_error(L, "realloc failed", ENOMEM);
    return NULL;
  } else {
    return p;
  }
}

static inline int tk_lua_absindex (lua_State *L, int i)
{
  if (i < 0 && i > LUA_REGISTRYINDEX)
    i += lua_gettop(L) + 1;
  return i;
}

static inline void tk_lua_register (lua_State *L, luaL_Reg *regs, int nup)
{
  while (true) {
    if ((*regs).name == NULL)
      break;
    for (int i = 0; i < nup; i ++)
      lua_pushvalue(L, -nup); // t upsa upsb
    lua_pushcclosure(L, (*regs).func, nup); // t upsa fn
    lua_setfield(L, -nup - 2, (*regs).name); // t
    regs ++;
  }
  lua_pop(L, nup);
}

static tk_compressor_t *peek_compressor (lua_State *L, int i)
{
  return (tk_compressor_t *) luaL_checkudata(L, i, MT_COMPRESSOR);
}

static inline void tk_compressor_shrink (tk_compressor_t *C)
{
  free(C->maxmis); C->maxmis = NULL;
  free(C->mis); C->mis = NULL;
  free(C->sums); C->sums = NULL;
  free(C->px); C->px = NULL;
  free(C->entropy_x); C->entropy_x = NULL;
  free(C->tmp); C->tmp = NULL;
  if (C->cnt) { kh_destroy(mtx, C->cnt); C->cnt = NULL; }
  free(C->smis); C->smis = NULL;
  free(C->tcs); C->tcs = NULL;
}

static int tk_compressor_gc (lua_State *L)
{
  lua_settop(L, 1);
  tk_compressor_t *C = peek_compressor(L, 1);
  if (C->destroyed)
    return 1;
  C->destroyed = true;
  tk_compressor_shrink(C);
  free(C->alpha); C->alpha = NULL;
  kh_destroy(mtx, C->lm); C->lm = NULL;
  free(C->log_py); C->log_py = NULL;
  free(C->log_pyx_unnorm); C->log_pyx_unnorm = NULL;
  free(C->pyx); C->pyx = NULL;
  free(C->sumsa); C->sumsa = NULL;
  free(C->samples); C->samples = NULL;
  free(C->visibles); C->visibles = NULL;
  pthread_mutex_lock(&C->mutex);
  C->stage = TK_CMP_DONE;
  pthread_cond_broadcast(&C->cond_stage);
  pthread_mutex_unlock(&C->mutex);
  // TODO: What is the right way to deal with potential thread errors (or other
  // errors, for that matter) during the finalizer?
  if (C->created_threads)
    for (unsigned int i = 0; i < C->n_threads; i ++)
      pthread_join(C->threads[i], NULL);
    // if (pthread_join(C->threads[i], NULL) != 0)
    //   tk_error(L, "pthread_join", errno);
  pthread_mutex_destroy(&C->mutex);
  pthread_cond_destroy(&C->cond_stage);
  pthread_cond_destroy(&C->cond_done);
  free(C->threads); C->threads = NULL;
  free(C->thread_data); C->thread_data = NULL;
  return 0;
}

static inline int tk_compressor_destroy (lua_State *L)
{
  lua_settop(L, 0);
  lua_pushvalue(L, lua_upvalueindex(1));
  return tk_compressor_gc(L);
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

static inline void seed_rand ()
{
  mcg_state = (uint64_t) pthread_self() ^ (uint64_t) time(NULL);
}

static inline int tk_compressor_compress (lua_State *);

static inline void tk_compressor_signal (
  tk_compressor_stage_t stage,
  tk_compressor_stage_t *stagep,
  pthread_mutex_t *mutex,
  pthread_cond_t *cond_stage,
  pthread_cond_t *cond_done,
  unsigned int *n_threads_done,
  unsigned int n_threads
) {
  pthread_mutex_lock(mutex);
  (*stagep) = stage;
  (*n_threads_done) = 0;
  pthread_cond_broadcast(cond_stage);
  pthread_mutex_unlock(mutex);
  pthread_mutex_lock(mutex);
  while ((*n_threads_done) < n_threads)
    pthread_cond_wait(cond_done, mutex);
  pthread_cond_broadcast(cond_stage);
  pthread_mutex_unlock(mutex);
}

static inline void tk_compressor_marginals_thread (
  tk_compressor_t *C,
  uint64_t cardinality,
  uint64_t *restrict samples,
  unsigned int *restrict visibles,
  double *restrict log_py,
  double *restrict pyx,
  double *restrict tmp,
  double *restrict mis,
  double *restrict smis,
  double *restrict px,
  double *restrict entropy_x,
  double *restrict tcs,
  double tmin,
  double ttc,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  double *restrict tmp00 = tmp + 0 * n_visible;
  double *restrict tmp01 = tmp + 1 * n_visible;
  double *restrict tmp10 = tmp + 2 * n_visible;
  double *restrict tmp11 = tmp + 3 * n_visible;
  double *restrict lpy0 = log_py + 0 * n_hidden;
  double *restrict lpy1 = log_py + 1 * n_hidden;
  double *restrict smis0 = smis + 0 * n_hidden * n_visible;
  double *restrict smis1 = smis + 1 * n_hidden * n_visible;
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    double *restrict pyx0 = pyx + 0 * n_hidden * n_samples + h * n_samples;
    double *restrict pyx1 = pyx + 1 * n_hidden * n_samples + h * n_samples;
    double counts_0 = 0.001;
    double counts_1 = 0.001;
    for (unsigned int s = 0; s < n_samples; s ++) {
      counts_0 += pyx0[s];
      counts_1 += pyx1[s];
    }
    double sum_counts = counts_0 + counts_1;
    lpy0[h] = log(counts_0) - log(sum_counts);
    lpy1[h] = log(counts_1) - log(sum_counts);
  }
  pc_clear(C);
  // for (unsigned int h = hfirst; h <= hlast; h ++) {
  //   double sum_py0 = 0.0;
  //   double sum_py1 = 0.0;
  //   double *restrict hpy0s = pyx + 0 * n_hidden * n_samples + h * n_samples;
  //   double *restrict hpy1s = pyx + 1 * n_hidden * n_samples + h * n_samples;
  //   for (unsigned int s = 0; s < n_samples; s ++)
  //     sum_py0 += hpy0s[s];
  //   for (unsigned int s = 0; s < n_samples; s ++)
  //     sum_py1 += hpy1s[s];
  //   for (unsigned int v = 0; v < n_visible; v ++) {
  //     pc_add(C, 0, 0, h, v, sum_py0);
  //     pc_add(C, 0, 1, h, v, sum_py1);
  //   }
  // }
    fprintf(stderr, "> a cnt:  %u\n", kh_size(C->cnt));
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    for (uint64_t c = 0; c < cardinality; c ++) {
      uint64_t s = samples[c];
      unsigned int v = visibles[c];
      double py0 = pyx[0 * n_hidden * n_samples + h * n_samples + s];
      double py1 = pyx[1 * n_hidden * n_samples + h * n_samples + s];
      pc_add(C, 1, 0, h, v, py0);
      pc_add(C, 1, 1, h, v, py1);
      pc_add(C, 0, 0, h, v, -py0);
      pc_add(C, 0, 1, h, v, -py1);
    }
  }
    fprintf(stderr, "> b cnt:  %u\n", kh_size(C->cnt));
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    for (unsigned int v = 0; v < n_visible; v ++) {
      double pc00 = pc_get(C, 0, 0, h, v) + 0.001;
      double pc01 = pc_get(C, 0, 1, h, v) + 0.001;
      double pc10 = pc_get(C, 1, 0, h, v) + 0.001;
      double pc11 = pc_get(C, 1, 1, h, v) + 0.001;
      double log_total0 = log(pc00 + pc01);
      double log_total1 = log(pc10 + pc11);
      lm_set(C, 0, 0, h, v, log(pc00) - log_total0);
      lm_set(C, 0, 1, h, v, log(pc01) - log_total0);
      lm_set(C, 1, 0, h, v, log(pc10) - log_total1);
      lm_set(C, 1, 1, h, v, log(pc11) - log_total1);
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    double lpy0v = lpy0[h];
    double lpy1v = lpy1[h];
    for (unsigned int v = 0; v < n_visible; v ++) {
      lm_add(C, 0, 0, h, v, -lpy0v);
      lm_add(C, 0, 1, h, v, -lpy1v);
      lm_add(C, 1, 0, h, v, -lpy0v);
      lm_add(C, 1, 1, h, v, -lpy1v);
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, non-affine base
    double lpy0v = lpy0[h];
    double lpy1v = lpy1[h];
    for (unsigned int v = 0; v < n_visible; v ++) {
      pc_add(C, 0, 0, h, v, exp(lm_get(C, 0, 0, h, v) + lpy0v));
      pc_add(C, 0, 1, h, v, exp(lm_get(C, 0, 1, h, v) + lpy1v));
      pc_add(C, 1, 0, h, v, exp(lm_get(C, 1, 0, h, v) + lpy0v));
      pc_add(C, 1, 1, h, v, exp(lm_get(C, 1, 1, h, v) + lpy1v));
    }
  }
    fprintf(stderr, "> c cnt:  %u\n", kh_size(C->cnt));
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    for (unsigned int v = 0; v < n_visible; v ++) {
      double pc00 = pc_get(C, 0, 0, h, v) + 0.001;
      double pc01 = pc_get(C, 0, 1, h, v) + 0.001;
      double pc10 = pc_get(C, 1, 0, h, v) + 0.001;
      double pc11 = pc_get(C, 1, 1, h, v) + 0.001;
      smis0[h * n_visible + v] = pc00 * lm_get(C, 0, 0, h, v) + pc01 * lm_get(C, 0, 1, h, v);
      smis1[h * n_visible + v] = pc10 * lm_get(C, 1, 0, h, v) + pc11 * lm_get(C, 1, 1, h, v);
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, dominance boundary
    double *restrict smis0a = smis0 + h * n_visible;
    double *restrict smis1a = smis1 + h * n_visible;
    double *restrict misa = mis + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      double weighted_sum = smis0a[v] * (1 - px[v]) + smis1a[v] * (px[v]);
      misa[v] = weighted_sum;
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, non-affine base
    double *restrict mish = mis + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      mish[v] /= entropy_x[v];
  }
  for (unsigned int h = hfirst; h <= hlast; h ++)
    tcs[h] = fabs(tcs[h]) * ttc + tmin;
}

static inline void tk_compressor_maxmis_thread (
  double *restrict mis,
  double *restrict maxmis,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int vfirst,
  unsigned int vlast
) {
  for (unsigned int v = vfirst; v <= vlast; v ++) { // not vectorized, unsupported outer form
    double max_val = 0.0;
    for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, costings not worth while
      double candidate = mis[h * n_visible + v];
      if (candidate > max_val)
        max_val = candidate;
    }
    maxmis[v] = max_val;
  }
}

static inline void tk_compressor_alpha_thread (
  tk_compressor_t *C,
  double *restrict alpha,
  double *restrict sumsa,
  double *restrict tcs,
  double *restrict mis,
  double *restrict maxmis,
  double lam,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  double *restrict sumsa0 = sumsa + 0 * n_hidden;
  double *restrict sumsa1 = sumsa + 1 * n_hidden;
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, non-affine
    double *restrict alphah = alpha + h * n_visible;
    double *restrict mish = mis + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      alphah[v] = (1.0 - lam) * alphah[v] + lam * exp(tcs[h] * (mish[v] - maxmis[v]));
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, unsupported outer form
    double s0 = 0.0, s1 = 0.0;
    double *restrict aph = alpha + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      s0 += aph[v] * lm_get(C, 0, 0, h, v);
      s1 += aph[v] * lm_get(C, 0, 1, h, v);
    }
    sumsa0[h] = s0;
    sumsa1[h] = s1;
  }
}

static inline void tk_compressor_latent_sumsa_thread (
  double *restrict sums,
  double *restrict sumsa,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  double *restrict sums0 = sums + 0 * n_hidden * n_samples;
  double *restrict sums1 = sums + 1 * n_hidden * n_samples;
  double *restrict sumsa0 = sumsa + 0 * n_hidden;
  double *restrict sumsa1 = sumsa + 1 * n_hidden;
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, not affine
    double s0 = sumsa0[h];
    double s1 = sumsa1[h];
    double *restrict sums0a = sums0 + h * n_samples;
    double *restrict sums1a = sums1 + h * n_samples;
    for (unsigned int i = 0; i < n_samples; i ++) {
      sums0a[i] = s0;
      sums1a[i] = s1;
    }
  }
}

static inline void tk_compressor_latent_sums_thread (
  tk_compressor_t *C,
  uint64_t *restrict samples,
  unsigned int *restrict visibles,
  double *restrict alpha,
  double *restrict sums,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int bfirst,
  unsigned int blast
) {
  double *restrict sums0 = sums + 0 * n_hidden * n_samples;
  double *restrict sums1 = sums + 1 * n_hidden * n_samples;
  for (unsigned int b = bfirst; b <= blast; b ++) {
    uint64_t s = samples[b];
    unsigned int v = visibles[b];
    for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, gather/scatter
      double *restrict sums0h = sums0 + h * n_samples;
      double *restrict sums1h = sums1 + h * n_samples;
      double *restrict aph = alpha + h * n_visible;
      sums0h[s] = sums0h[s] - aph[v] * lm_get(C, 0, 0, h, v) + aph[v] * lm_get(C, 1, 0, h, v); // gather/scatter
      sums1h[s] = sums1h[s] - aph[v] * lm_get(C, 0, 1, h, v) + aph[v] * lm_get(C, 1, 1, h, v); // gather/scatter
    }
  }
}

static inline void tk_compressor_latent_py_thread (
  double *restrict log_py,
  double *restrict log_pyx_unnorm,
  double *restrict sums,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  double *restrict sums0 = sums + 0 * n_hidden * n_samples;
  double *restrict sums1 = sums + 1 * n_hidden * n_samples;
  double *restrict lpy0 = log_py + 0 * n_hidden;
  double *restrict lpy1 = log_py + 1 * n_hidden;
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, not affine
    double lpy0v = lpy0[h];
    double lpy1v = lpy1[h];
    double *restrict sums0h = sums0 + h * n_samples;
    double *restrict sums1h = sums1 + h * n_samples;
    double *restrict lpyx0 = log_pyx_unnorm + 0 * n_hidden * n_samples + h * n_samples;
    double *restrict lpyx1 = log_pyx_unnorm + 1 * n_hidden * n_samples + h * n_samples;
    for (unsigned int i = 0; i < n_samples; i ++) {
      lpyx0[i] = sums0h[i] + lpy0v;
      lpyx1[i] = sums1h[i] + lpy1v;
    }
  }
}

static inline void tk_compressor_latent_norm_thread (
  double *restrict log_z, // mis
  double *restrict pyx,
  double *restrict log_pyx_unnorm,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  double *restrict pyx0 = pyx + 0 * n_hidden * n_samples;
  double *restrict pyx1 = pyx + 1 * n_hidden * n_samples;
  double *restrict lpyx0 = log_pyx_unnorm + 0 * n_hidden * n_samples;
  double *restrict lpyx1 = log_pyx_unnorm + 1 * n_hidden * n_samples;
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++) {
    double a = lpyx0[i];
    double b = lpyx1[i];
    double max_ab  = (a > b) ? a : b;
    double sum_exp = exp(a - max_ab) + exp(b - max_ab);
    log_z[i] = max_ab + log(sum_exp);
  }
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++)
    pyx0[i] = exp(lpyx0[i] - log_z[i]);
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++)
    pyx1[i] = exp(lpyx1[i] - log_z[i]);
}

static inline void tk_compressor_update_tc_thread (
  double *restrict log_z, // mis
  double *restrict tcs,
  unsigned int n_samples,
  unsigned int hfirst,
  unsigned int hlast
) {
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, unsupported outer form
    const double *restrict lz = log_z + h * n_samples;
    double sum = 0.0;
    for (unsigned int s = 0; s < n_samples; s ++)
      sum += lz[s];
    double tc = sum / (double) n_samples;
    tcs[h] = tc;
  }
}

static inline void tk_compressor_update_last_tc (
  double *restrict tcs,
  double *last_tc,
  unsigned int n_hidden
) {
  double sum0 = 0.0;
  for (unsigned int h = 0; h < n_hidden; h ++)
    sum0 += tcs[h];
  *last_tc = sum0;
}

typedef struct {
  uint64_t *samples;
  unsigned int *visibles;
  uint64_t n;
  unsigned int n_visible;
} tk_compressor_setup_bits_t;

static bool tk_compressor_setup_bits_iter (uint64_t val, void *statepv)
{
  tk_compressor_setup_bits_t *statep =
    (tk_compressor_setup_bits_t *) statepv;
  statep->samples[statep->n] = val / statep->n_visible;
  statep->visibles[statep->n] = val % statep->n_visible;
  statep->n ++;
  return true;
}

static inline void tk_compressor_setup_bits (
  roaring64_bitmap_t *bm,
  uint64_t *samples,
  unsigned int *visibles,
  unsigned int n_visible
) {
  tk_compressor_setup_bits_t state;
  state.samples = samples;
  state.visibles = visibles;
  state.n = 0;
  state.n_visible = n_visible;
  roaring64_bitmap_iterate(bm, tk_compressor_setup_bits_iter, &state);
}

// Does this make sense to parallize? I don't think so...
static inline void tk_compressor_data_stats (
  uint64_t cardinality,
  unsigned int *restrict visibles,
  double *restrict px,
  double *restrict entropy_x,
  unsigned int n_samples,
  unsigned int n_visible
) {
  for (unsigned int v = 0; v < n_visible; v ++)
    px[v] = 0;
  for (uint64_t c = 0; c < cardinality; c ++)
    px[visibles[c]] ++;
  for (unsigned int v = 0; v < n_visible; v ++)
    px[v] /= (double) n_samples;
  for (unsigned int v = 0; v < n_visible; v ++) {
    double entropy = 0;
    entropy -= px[v] * log(px[v]);
    entropy -= (1 - px[v]) * log(1 - px[v]);
    entropy_x[v] = entropy > 0 ? entropy : 0.001;
  }
}

static inline void tk_compressor_init_alpha_thread (
  double *alpha,
  unsigned int n_visible,
  unsigned int hfirst,
  unsigned int hlast
) {
  for (unsigned int i = hfirst * n_visible; i < (hlast + 1) * n_visible; i++)
    alpha[i] = 0.5 + 0.5 * fast_drand();
}

// Doesn't really need to be threaded..'
// Consider combining with one of the other thread inits
static inline void tk_compressor_init_tcs_thread (
  double *tcs,
  unsigned int hfirst,
  unsigned int hlast
) {
  for (unsigned int i = hfirst; i <= hlast; i ++)
    tcs[i] = 0.0;
}

static inline void tk_compressor_init_log_pyx_unnorm_thread (
  double *log_pyx_unnorm,
  unsigned int n_samples,
  unsigned int hfirst,
  unsigned int hlast
) {
  double log_dim_hidden = -log(2);
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++)
    log_pyx_unnorm[i] = log_dim_hidden * (0.5 + fast_drand());
}

static void *tk_compressor_worker (void *datap)
{
  seed_rand();
  tk_compressor_thread_data_t *data =
    (tk_compressor_thread_data_t *) datap;
  while (1) {
    pthread_mutex_lock(&data->C->mutex);
    while (data->stage == data->C->stage)
      pthread_cond_wait(&data->C->cond_stage, &data->C->mutex);
    data->stage = data->C->stage;
    pthread_mutex_unlock(&data->C->mutex);
    if (data->stage == TK_CMP_DONE)
      break;
    switch (data->stage) {
      case TK_CMP_INIT_ALPHA:
        tk_compressor_init_alpha_thread(
          data->C->alpha,
          data->C->n_visible,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_INIT_TCS:
        tk_compressor_init_tcs_thread(
          data->C->tcs,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_INIT_PYX_UNNORM:
        tk_compressor_init_log_pyx_unnorm_thread(
          data->C->log_pyx_unnorm,
          data->n_samples,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_MARGINALS:
        tk_compressor_marginals_thread(
          data->C,
          data->cardinality,
          data->C->samples,
          data->C->visibles,
          data->C->log_py,
          data->C->pyx,
          data->C->tmp,
          data->C->mis,
          data->C->smis,
          data->C->px,
          data->C->entropy_x,
          data->C->tcs,
          data->C->tmin,
          data->C->ttc,
          data->n_samples,
          data->C->n_visible,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_MAXMIS:
        tk_compressor_maxmis_thread(
          data->C->mis,
          data->C->maxmis,
          data->C->n_visible,
          data->C->n_hidden,
          data->vfirst,
          data->vlast);
        break;
      case TK_CMP_ALPHA:
        tk_compressor_alpha_thread(
          data->C,
          data->C->alpha,
          data->C->sumsa,
          data->C->tcs,
          data->C->mis,
          data->C->maxmis,
          data->C->lam,
          data->C->n_visible,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_LATENT_SUMSA:
        tk_compressor_latent_sumsa_thread(
          data->C->sums,
          data->C->sumsa,
          data->n_samples,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_LATENT_SUMS:
        tk_compressor_latent_sums_thread(
          data->C,
          data->C->samples,
          data->C->visibles,
          data->C->alpha,
          data->C->sums,
          data->n_samples,
          data->C->n_visible,
          data->C->n_hidden,
          data->bfirst,
          data->blast);
        break;
      case TK_CMP_LATENT_PY:
        tk_compressor_latent_py_thread(
          data->C->log_py,
          data->C->log_pyx_unnorm,
          data->C->sums,
          data->n_samples,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_LATENT_NORM:
        tk_compressor_latent_norm_thread(
          data->C->mis,
          data->C->pyx,
          data->C->log_pyx_unnorm,
          data->n_samples,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_UPDATE_TC:
        tk_compressor_update_tc_thread(
          data->C->mis,
          data->C->tcs,
          data->n_samples,
          data->hfirst,
          data->hlast);
        break;
      default:
        assert(false);
    }
    pthread_mutex_lock(&data->C->mutex);
    data->C->n_threads_done ++;
    if (data->C->n_threads_done == data->C->n_threads)
      pthread_cond_signal(&data->C->cond_done);
    pthread_mutex_unlock(&data->C->mutex);
  }
  return NULL;
}

static inline void tk_compressor_setup_threads (
  lua_State *L,
  tk_compressor_t *C,
  unsigned int cardinality,
  unsigned int n_samples
) {
  unsigned int hslice = C->n_hidden / C->n_threads;
  unsigned int hremaining = C->n_hidden % C->n_threads;
  unsigned int hfirst = 0;
  unsigned int vslice = C->n_visible / C->n_threads;
  unsigned int vremaining = C->n_visible % C->n_threads;
  unsigned int vfirst = 0;
  unsigned int bslice = cardinality / C->n_threads;
  unsigned int bremaining = cardinality % C->n_threads;
  unsigned int bfirst = 0;
  for (unsigned int i = 0; i < C->n_threads; i ++) {
    C->thread_data[i].C = C;
    C->thread_data[i].cardinality = cardinality;
    C->thread_data[i].n_samples = n_samples;
    C->thread_data[i].stage = TK_CMP_INIT;
    C->thread_data[i].hfirst = hfirst;
    C->thread_data[i].hlast = hfirst + hslice - 1;
    if (hremaining) {
      C->thread_data[i].hlast ++;
      hremaining --;
    }
    hfirst = C->thread_data[i].hlast + 1;
    C->thread_data[i].vfirst = vfirst;
    C->thread_data[i].vlast = vfirst + vslice - 1;
    if (vremaining) {
      C->thread_data[i].vlast ++;
      vremaining --;
    }
    vfirst = C->thread_data[i].vlast + 1;
    C->thread_data[i].bfirst = bfirst;
    C->thread_data[i].blast = bfirst + bslice - 1;
    if (bremaining) {
      C->thread_data[i].blast ++;
      bremaining --;
    }
    while (C->thread_data[i].blast + 1 < cardinality &&
      C->samples[C->thread_data[i].blast + 1] == C->samples[C->thread_data[i].blast])
      C->thread_data[i].blast ++;
    if (C->thread_data[i].blast >= cardinality)
      C->thread_data[i].blast = cardinality - 1;
    bfirst = C->thread_data[i].blast + 1;
    // TODO: ensure everything gets freed on error (should be in compressor gc)
    if (!C->created_threads && pthread_create(&C->threads[i], NULL, tk_compressor_worker, &C->thread_data[i]) != 0)
      tk_error(L, "pthread_create", errno);
  }
  C->created_threads = true;
}

static inline int tk_compressor_compress (lua_State *L)
{
  tk_compressor_t *C = peek_compressor(L, lua_upvalueindex(1));
  roaring64_bitmap_t *bm = peek_bitmap(L, 1);
  uint64_t cardinality = roaring64_bitmap_get_cardinality(bm);
  // TODO: Expose shrink via the api, and only realloc if new size is larger than old
  C->samples = tk_realloc(L, C->samples, cardinality * sizeof(uint64_t));
  C->visibles = tk_realloc(L, C->visibles, cardinality * sizeof(unsigned int));
  tk_compressor_setup_bits(bm, C->samples, C->visibles, C->n_visible);
  unsigned int n_samples = tk_lua_optunsigned(L, 2, 1);
  tk_compressor_setup_threads(L, C, cardinality, n_samples);
  C->mis = tk_realloc(L, C->mis, C->n_hidden * n_samples * sizeof(double));
  C->pyx = tk_realloc(L, C->pyx, 2 * C->n_hidden * n_samples * sizeof(double));
  C->log_pyx_unnorm = tk_realloc(L, C->log_pyx_unnorm, 2 * C->n_hidden * n_samples * sizeof(double));
  unsigned int len_sums = (2 * C->n_hidden * (n_samples > C->n_visible ? n_samples : C->n_visible));
  C->sums = tk_realloc(L, C->sums, len_sums * sizeof(double));
  tk_compressor_signal(
    TK_CMP_LATENT_SUMSA,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  tk_compressor_signal(
    TK_CMP_LATENT_SUMS,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  tk_compressor_signal(
    TK_CMP_LATENT_PY,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  tk_compressor_signal(
    TK_CMP_LATENT_NORM,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  roaring64_bitmap_t *bm0 = roaring64_bitmap_create();
  if (bm0 == NULL)
    luaL_error(L, "memory error creating bitmap");
  roaring64_bitmap_t **bm0p = (roaring64_bitmap_t **)
    lua_newuserdata(L, sizeof(roaring64_bitmap_t *));
  *bm0p = bm0;
  luaL_getmetatable(L, MT_BITMAP);
  lua_setmetatable(L, -2);
  // TODO: Can this be parallelized? Is it worth it?
  for (unsigned int h = 0; h < C->n_hidden; h ++) {
    for (unsigned int s = 0; s < n_samples; s ++) {
      double py0 = C->pyx[0 * C->n_hidden * n_samples + h * n_samples + s];
      double py1 = C->pyx[1 * C->n_hidden * n_samples + h * n_samples + s];
      if (py1 > py0)
        roaring64_bitmap_add(bm0, s * C->n_hidden + h);
    }
  }
  return 1;
}

static inline void _tk_compressor_train (
  lua_State *L,
  tk_compressor_t *C,
  roaring64_bitmap_t *bm,
  unsigned int n_samples,
  unsigned int max_iter,
  int i_each
) {
  C->pyx = tk_malloc(L, 2 * C->n_hidden * n_samples * sizeof(double));
  C->log_pyx_unnorm = tk_malloc(L, 2 * C->n_hidden * n_samples * sizeof(double));
  unsigned int len_mis = C->n_hidden * C->n_visible;
  if (len_mis < (C->n_hidden * n_samples))
    len_mis = C->n_hidden * n_samples;
  C->mis = tk_malloc(L, len_mis * sizeof(double));
  C->sums = tk_malloc(L, 2 * n_samples * C->n_hidden * sizeof(double));
  uint64_t cardinality = roaring64_bitmap_get_cardinality(bm);
  C->samples = tk_malloc(L, cardinality * sizeof(uint64_t));
  C->visibles = tk_malloc(L, cardinality * sizeof(unsigned int));
  tk_compressor_setup_bits(bm, C->samples, C->visibles, C->n_visible);
  tk_compressor_setup_threads(L, C, cardinality, n_samples);
  tk_compressor_data_stats(
    cardinality,
    C->visibles,
    C->px,
    C->entropy_x,
    n_samples,
    C->n_visible);
  tk_compressor_signal(
    TK_CMP_INIT_TCS,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  tk_compressor_signal(
    TK_CMP_INIT_ALPHA,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  tk_compressor_signal(
    TK_CMP_INIT_PYX_UNNORM,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  tk_compressor_signal(
    TK_CMP_LATENT_NORM,
    &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
    &C->n_threads_done, C->n_threads);
  for (unsigned int i = 0; i < max_iter; i ++) {
    tk_compressor_signal(
      TK_CMP_MARGINALS,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_MAXMIS,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_ALPHA,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_LATENT_SUMSA,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_LATENT_SUMS,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_LATENT_PY,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_LATENT_NORM,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_signal(
      TK_CMP_UPDATE_TC,
      &C->stage, &C->mutex, &C->cond_stage, &C->cond_done,
      &C->n_threads_done, C->n_threads);
    tk_compressor_update_last_tc(
      C->tcs,
      &C->last_tc,
      C->n_hidden);
    if (i_each > -1) {
      lua_pushvalue(L, i_each);
      lua_pushinteger(L, i + 1);
      lua_pushnumber(L, C->last_tc);
      lua_call(L, 2, 1);
      if (lua_type(L, -1) == LUA_TBOOLEAN && lua_toboolean(L, -1) == 0)
        break;
      lua_pop(L, 1);
    }
  }
  tk_compressor_shrink(C);
  C->trained = true;
}

static inline void tk_compressor_init (
  lua_State *L,
  tk_compressor_t *C,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int n_threads
) {
  memset(C, 0, sizeof(tk_compressor_t));
  C->n_visible = n_visible;
  C->n_hidden = n_hidden;
  C->lam = 0.3;
  C->tmin = 1.0;
  C->ttc = 500.0;
  C->tcs = tk_malloc(L, C->n_hidden * sizeof(double));
  C->log_py = tk_malloc(L, 2 * C->n_hidden * sizeof(double));
  C->lm = kh_init(mtx);
  C->cnt = kh_init(mtx);
  C->alpha = tk_malloc(L, C->n_hidden * C->n_visible * sizeof(double));
  C->tmp = tk_malloc(L, 2 * 2 * C->n_visible * sizeof(double));
  C->smis = tk_malloc(L, 2 * C->n_hidden * C->n_visible * sizeof(double));
  C->px = tk_malloc(L, C->n_visible * sizeof(double));
  C->entropy_x = tk_malloc(L, C->n_visible * sizeof(double));
  C->maxmis = tk_malloc(L, C->n_visible * sizeof(double));
  C->sumsa = tk_malloc(L, 2 * C->n_hidden * sizeof(double));
  C->n_threads = 1; //n_threads;
  C->n_threads_done = 0;
  C->stage = TK_CMP_INIT;
  C->threads = tk_malloc(L, C->n_threads * sizeof(pthread_t));
  C->thread_data = tk_malloc(L, C->n_threads * sizeof(tk_compressor_thread_data_t));
  // TODO: check errors
  pthread_mutex_init(&C->mutex, NULL);
  pthread_cond_init(&C->cond_stage, NULL);
  pthread_cond_init(&C->cond_done, NULL);
}

static inline int tk_compressor_train (lua_State *L) {
  tk_compressor_t *C = peek_compressor(L, lua_upvalueindex(1));
  if (C->trained)
    return tk_lua_error(L, "Already trained!\n");
  roaring64_bitmap_t *bm = *((roaring64_bitmap_t **)
    tk_lua_fcheckuserdata(L, 1, "corpus", MT_BITMAP));
  unsigned int n_samples = tk_lua_fcheckunsigned(L, 1, "samples");
  unsigned int max_iter = tk_lua_fcheckunsigned(L, 1, "iterations");
  int i_each = -1;
  if (tk_lua_ftype(L, 1, "each") != LUA_TNIL) {
    lua_getfield(L, 1, "each");
    i_each = tk_lua_absindex(L, -1);
  }
  _tk_compressor_train(L, C, bm, n_samples, max_iter, i_each); // c
  return 0;
}

static inline int tk_compressor_persist (lua_State *L)
{
  tk_compressor_t *C = peek_compressor(L, lua_upvalueindex(1));
  if (!C->trained)
    return tk_lua_error(L, "Can't persist an untrained model\n");
  lua_settop(L, 1);
  bool tostr = lua_type(L, 1) == LUA_TNIL;
  FILE *fh = tostr ? tk_lua_tmpfile(L) : tk_lua_fopen(L, luaL_checkstring(L, 1), "w");
  tk_lua_fwrite(L, &C->trained, sizeof(bool), 1, fh);
  tk_lua_fwrite(L, &C->n_visible, sizeof(unsigned int), 1, fh);
  tk_lua_fwrite(L, &C->n_hidden, sizeof(unsigned int), 1, fh);
  tk_lua_fwrite(L, &C->lam, sizeof(double), 1, fh);
  tk_lua_fwrite(L, &C->tmin, sizeof(double), 1, fh);
  tk_lua_fwrite(L, &C->ttc, sizeof(double), 1, fh);
  tk_lua_fwrite(L, C->alpha, sizeof(double), C->n_hidden * C->n_visible, fh);
  tk_lua_fwrite(L, C->log_py, sizeof(double), 2 * C->n_hidden, fh);
  #warning "Write lm to disk"
  tk_lua_fwrite(L, C->sumsa, sizeof(double), 2 * C->n_hidden, fh);
  if (!tostr) {
    tk_lua_fclose(L, fh);
    return 0;
  } else {
    size_t len;
    char *data = tk_lua_fslurp(L, fh, &len);
    if (data) {
      lua_pushlstring(L, data, len);
      free(data);
      tk_lua_fclose(L, fh);
      return 1;
    } else {
      tk_lua_fclose(L, fh);
      return 0;
    }
  }
}

static luaL_Reg mt_fns[] =
{
  { "compress", tk_compressor_compress },
  { "persist", tk_compressor_persist },
  { "train", tk_compressor_train },
  { "destroy", tk_compressor_destroy },
  { NULL, NULL }
};

static inline int tk_compressor_load (lua_State *L)
{
  lua_settop(L, 3); // fp ts
  tk_compressor_t *C = (tk_compressor_t *)
    lua_newuserdata(L, sizeof(tk_compressor_t)); // tp ts c
  memset(C, 0, sizeof(tk_compressor_t));
  luaL_getmetatable(L, MT_COMPRESSOR); // fp ts c mt
  lua_setmetatable(L, -2); // fp ts c
  lua_newtable(L); // fp ts c t
  lua_pushvalue(L, -2); // fp ts c t c
  tk_lua_register(L, mt_fns, 1); // fp ts c t
  lua_remove(L, -2); // fp ts t
  unsigned int n_threads;
  if (lua_type(L, 2) != LUA_TNIL) {
    // TODO: allow passing 0 to run everything on the main thread.
    n_threads = tk_lua_checkunsigned(L, 2);
    if (!n_threads)
      return tk_lua_error(L, "threads must be at least 1\n");
  } else {
    long ts = sysconf(_SC_NPROCESSORS_ONLN) - 1;
    if (ts <= 0)
      return tk_error(L, "sysconf", errno);
    lua_pushinteger(L, ts);
    n_threads = tk_lua_checkunsigned(L, -1);
    lua_pop(L, 1);
  }
  size_t len;
  const char *data = luaL_checklstring(L, 1, &len);
  bool isstr = lua_type(L, 3) == LUA_TBOOLEAN && lua_toboolean(L, 3);
  FILE *fh = isstr ? tk_lua_fmemopen(L, (char *) data, len, "r") : tk_lua_fopen(L, data, "r");
  tk_lua_fread(L, &C->trained, sizeof(bool), 1, fh);
  tk_lua_fread(L, &C->n_visible, sizeof(unsigned int), 1, fh);
  tk_lua_fread(L, &C->n_hidden, sizeof(unsigned int), 1, fh);
  tk_lua_fread(L, &C->lam, sizeof(double), 1, fh);
  tk_lua_fread(L, &C->tmin, sizeof(double), 1, fh);
  tk_lua_fread(L, &C->ttc, sizeof(double), 1, fh);
  C->alpha = tk_malloc(L, C->n_hidden * C->n_visible * sizeof(double));
  tk_lua_fread(L, C->alpha, sizeof(double), C->n_hidden * C->n_visible, fh);
  C->log_py = tk_malloc(L, 2 * C->n_hidden * sizeof(double));
  tk_lua_fread(L, C->log_py, sizeof(double), 2 * C->n_hidden, fh);
  C->lm = kh_init(mtx);
  #warning "Read btree lm from disk"
  C->tmp = tk_malloc(L, 2 * 2 * C->n_visible * sizeof(double));
  C->sumsa = tk_malloc(L, 2 * C->n_hidden * sizeof(double));
  tk_lua_fread(L, C->sumsa, sizeof(double), 2 * C->n_hidden, fh);
  C->n_threads = n_threads;
  C->n_threads_done = 0;
  C->stage = TK_CMP_INIT;
  C->threads = tk_malloc(L, C->n_threads * sizeof(pthread_t));
  C->thread_data = tk_malloc(L, C->n_threads * sizeof(tk_compressor_thread_data_t));
  tk_lua_fclose(L, fh);
  pthread_mutex_init(&C->mutex, NULL);
  pthread_cond_init(&C->cond_stage, NULL);
  pthread_cond_init(&C->cond_done, NULL);
  return 1;
}

static inline int tk_compressor_create (lua_State *L)
{
  lua_settop(L, 1);
  unsigned int n_visible = tk_lua_fcheckunsigned(L, 1, "visible");
  unsigned int n_hidden = tk_lua_fcheckunsigned(L, 1, "hidden");
  unsigned int n_threads;
  if (tk_lua_ftype(L, 1, "threads") != LUA_TNIL) {
    // TODO: allow passing 0 to run everything on the main thread.
    n_threads = tk_lua_fcheckunsigned(L, 1, "threads");
    if (!n_threads)
      return tk_lua_error(L, "threads must be at least 1\n");
  } else {
    long ts = sysconf(_SC_NPROCESSORS_ONLN) - 1;
    if (ts <= 0)
      return tk_error(L, "sysconf", errno);
    lua_pushinteger(L, ts);
    n_threads = tk_lua_checkunsigned(L, -1);
    lua_pop(L, 1);
  }
  tk_compressor_t *C = (tk_compressor_t *)
    lua_newuserdata(L, sizeof(tk_compressor_t)); // c
  luaL_getmetatable(L, MT_COMPRESSOR); // c mt
  lua_setmetatable(L, -2); // c
  tk_compressor_init(L, C, n_visible, n_hidden, n_threads); // c
  lua_newtable(L); // c t
  lua_pushvalue(L, -2); // c t c
  tk_lua_register(L, mt_fns, 1); // t
  return 1;
}

static luaL_Reg fns[] =
{
  { "create", tk_compressor_create },
  { "load", tk_compressor_load },
  { NULL, NULL }
};

int luaopen_santoku_bitmap_compressor (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, fns); // t
  luaL_newmetatable(L, MT_COMPRESSOR); // t mt
  lua_pushcfunction(L, tk_compressor_gc); // t mt fn
  lua_setfield(L, -2, "__gc"); // t mt
  lua_pop(L, 1); // t
  return 1;
}
