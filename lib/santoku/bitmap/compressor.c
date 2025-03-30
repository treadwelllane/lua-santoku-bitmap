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

#define MT_BITMAP "santoku_bitmap"
#define MT_COMPRESSOR "santoku_bitmap_compressor"

typedef enum {
  TK_CMP_INIT,
  TK_CMP_INIT_ALPHA,
  TK_CMP_INIT_TCS,
  TK_CMP_INIT_PYX_UNNORM,
  TK_CMP_MARGINALS,
  TK_CMP_MAXMIS,
  TK_CMP_ALPHA,
  TK_CMP_LATENT_BASELINE,
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
  float *alpha;
  float *log_marg;
  float *log_py;
  float *log_pyx_unnorm;
  float *maxmis;
  float *mis;
  float *sums;
  float *baseline;
  uint64_t *samples;
  unsigned int *visibles;
  float *px;
  float *entropy_x;
  float *pyx;
  float *counts;
  float last_tc;
  float *tcs;
  float lam;
  float tmin;
  float ttc;
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

static inline lua_Number tk_lua_foptnumber (lua_State *L, int i, char *field, float d)
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
  free(C->counts); C->counts = NULL;
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
  free(C->log_marg); C->log_marg = NULL;
  free(C->log_py); C->log_py = NULL;
  free(C->log_pyx_unnorm); C->log_pyx_unnorm = NULL;
  free(C->pyx); C->pyx = NULL;
  free(C->baseline); C->baseline = NULL;
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

static inline float fast_drand ()
{
  return ((float)fast_rand()) / ((float)UINT32_MAX);
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
  uint64_t cardinality,
  uint64_t *restrict samples,
  unsigned int *restrict visibles,
  float *restrict log_py,
  float *restrict pyx,
  float *restrict counts,
  float *restrict log_marg,
  float *restrict mis,
  float *restrict px,
  float *restrict entropy_x,
  float *restrict tcs,
  float tmin,
  float ttc,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  float *restrict lm00 = log_marg + 0 * n_hidden * n_visible;
  float *restrict lm01 = log_marg + 1 * n_hidden * n_visible;
  float *restrict lm10 = log_marg + 2 * n_hidden * n_visible;
  float *restrict lm11 = log_marg + 3 * n_hidden * n_visible;
  float *restrict pc00 = counts + 0 * n_hidden * n_visible;
  float *restrict pc01 = counts + 1 * n_hidden * n_visible;
  float *restrict pc10 = counts + 2 * n_hidden * n_visible;
  float *restrict pc11 = counts + 3 * n_hidden * n_visible;
  float *restrict tmp00 = lm00 + hfirst * n_visible;
  float *restrict tmp01 = lm01 + hfirst * n_visible;
  float *restrict tmp10 = lm10 + hfirst * n_visible;
  float *restrict tmp11 = lm11 + hfirst * n_visible;
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    float *restrict pyx0 = pyx + h * n_samples;
    float counts_0 = 0.001;
    float counts_1 = 0.001;
    for (unsigned int s = 0; s < n_samples; s ++) {
      counts_0 += pyx0[s];
      counts_1 += (1 - pyx0[s]);
    }
    float sum_counts = counts_0 + counts_1;
    log_py[h] = logf(counts_0) - logf(sum_counts);
  }
  for (unsigned int i = hfirst * n_visible; i < (hlast + 1) * n_visible; i++)
    counts[i] = 0;
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    float sum_py0 = 0.0;
    float *restrict hpy0s = pyx + h * n_samples;
    for (unsigned int s = 0; s < n_samples; s ++)
      sum_py0 += hpy0s[s];
    float sum_py1 = n_samples - sum_py0;
    float *restrict pc00a = pc00 + h * n_visible;
    float *restrict pc01a = pc01 + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      pc00a[v] += sum_py0;
      pc01a[v] += sum_py1;
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    for (unsigned int v = 0; v < n_visible; v ++) {
      tmp10[v] = 0;
      tmp11[v] = 0;
      tmp00[v] = 0;
      tmp01[v] = 0;
    }
    for (uint64_t c = 0; c < cardinality; c ++) { // not vectorzed, gather scatter
      uint64_t s = samples[c];
      unsigned int v = visibles[c];
      float py0 = pyx[h * n_samples + s];
      float py1 = 1 - py0;
      tmp10[v] += py0;
      tmp11[v] += py1;
      tmp00[v] -= py0;
      tmp01[v] -= py1;
    }
    float *restrict pc00a = pc00 + h * n_visible;
    float *restrict pc01a = pc01 + h * n_visible;
    float *restrict pc10a = pc10 + h * n_visible;
    float *restrict pc11a = pc11 + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      pc10a[v] += tmp10[v];
      pc11a[v] += tmp11[v];
      pc00a[v] += tmp00[v];
      pc01a[v] += tmp01[v];
    }
  }
  for (unsigned int i = hfirst * n_visible; i < (hlast + 1) * n_visible; i ++) {
    pc00[i] += 0.001;
    pc01[i] += 0.001;
    pc10[i] += 0.001;
    pc11[i] += 0.001;
  }
  for (unsigned int i = hfirst * n_visible; i < (hlast + 1) * n_visible; i ++) {
    float log_total0 = logf(pc00[i] + pc01[i]);
    lm00[i] = logf(pc00[i]) - log_total0;
    lm01[i] = logf(pc01[i]) - log_total0;
  }
  for (unsigned int i = hfirst * n_visible; i < (hlast + 1) * n_visible; i ++) {
    float log_total1 = logf(pc10[i] + pc11[i]);
    lm10[i] = logf(pc10[i]) - log_total1;
    lm11[i] = logf(pc11[i]) - log_total1;
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) {
    float *restrict lm00a = lm00 + h * n_visible;
    float *restrict lm01a = lm01 + h * n_visible;
    float *restrict lm10a = lm10 + h * n_visible;
    float *restrict lm11a = lm11 + h * n_visible;
    float lpy0v = log_py[h];
    float lpy1v = log1pf(-expf(lpy0v));
    for (unsigned int v = 0; v < n_visible; v ++) {
      lm00a[v] -= lpy0v;
      lm01a[v] -= lpy1v;
      lm10a[v] -= lpy0v;
      lm11a[v] -= lpy1v;
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, non-affine base
    float *restrict lm00a = lm00 + h * n_visible;
    float *restrict lm01a = lm01 + h * n_visible;
    float *restrict lm10a = lm10 + h * n_visible;
    float *restrict lm11a = lm11 + h * n_visible;
    float *restrict pc00a = pc00 + h * n_visible;
    float *restrict pc01a = pc01 + h * n_visible;
    float *restrict pc10a = pc10 + h * n_visible;
    float *restrict pc11a = pc11 + h * n_visible;
    float lpy0v = log_py[h];
    float lpy1v = log1pf(-expf(lpy0v));
    for (unsigned int v = 0; v < n_visible; v ++) {
      pc00a[v] = expf(lm00a[v] + lpy0v);
      pc01a[v] = expf(lm01a[v] + lpy1v);
      pc10a[v] = expf(lm10a[v] + lpy0v);
      pc11a[v] = expf(lm11a[v] + lpy1v);
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h++) {
    float *restrict pc00a = pc00 + h * n_visible;
    float *restrict pc01a = pc01 + h * n_visible;
    float *restrict pc10a = pc10 + h * n_visible;
    float *restrict pc11a = pc11 + h * n_visible;
    float *restrict lm00a = lm00 + h * n_visible;
    float *restrict lm10a = lm10 + h * n_visible;
    float *restrict lm01a = lm01 + h * n_visible;
    float *restrict lm11a = lm11 + h * n_visible;
    float *restrict misa   = mis   + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v++) {
      float group0 = pc00a[v] * lm00a[v] + pc01a[v] * lm01a[v];
      float group1 = pc10a[v] * lm10a[v] + pc11a[v] * lm11a[v];
      misa[v] = group0 * (1 - px[v]) + group1 * px[v];
    }
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, non-affine base
    float *restrict mish = mis + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      mish[v] /= entropy_x[v];
  }
  for (unsigned int h = hfirst; h <= hlast; h ++)
    tcs[h] = fabsf(tcs[h]) * ttc + tmin;
}

static inline void tk_compressor_maxmis_thread (
  float *restrict mis,
  float *restrict maxmis,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int vfirst,
  unsigned int vlast
) {
  for (unsigned int v = vfirst; v <= vlast; v ++) { // not vectorized, unsupported outer form
    float max_val = 0.0;
    for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, costings not worth while
      float candidate = mis[h * n_visible + v];
      if (candidate > max_val)
        max_val = candidate;
    }
    maxmis[v] = max_val;
  }
}

static inline void tk_compressor_alpha_thread (
  float *restrict alpha,
  float *restrict baseline,
  float *restrict log_marg,
  float *restrict tcs,
  float *restrict mis,
  float *restrict maxmis,
  float lam,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  float *restrict baseline0 = baseline + 0 * n_hidden;
  float *restrict baseline1 = baseline + 1 * n_hidden;
  float *restrict lm00 = log_marg + 0 * n_hidden * n_visible;
  float *restrict lm01 = log_marg + 1 * n_hidden * n_visible;
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, non-affine
    float *restrict alphah = alpha + h * n_visible;
    float *restrict mish = mis + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++)
      alphah[v] = (1.0 - lam) * alphah[v] + lam * expf(tcs[h] * (mish[v] - maxmis[v]));
  }
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, unsupported outer form
    float s0 = 0.0, s1 = 0.0;
    float *restrict lm00a = lm00 + h * n_visible;
    float *restrict lm01a = lm01 + h * n_visible;
    float *restrict aph = alpha + h * n_visible;
    for (unsigned int v = 0; v < n_visible; v ++) {
      s0 += aph[v] * lm00a[v];
      s1 += aph[v] * lm01a[v];
    }
    baseline0[h] = s0;
    baseline1[h] = s1;
  }
}

static inline void tk_compressor_latent_baseline_thread (
  float *restrict sums,
  float *restrict baseline,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  float *restrict sums0 = sums + 0 * n_hidden * n_samples;
  float *restrict sums1 = sums + 1 * n_hidden * n_samples;
  float *restrict baseline0 = baseline + 0 * n_hidden;
  float *restrict baseline1 = baseline + 1 * n_hidden;
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, not affine
    float s0 = baseline0[h];
    float s1 = baseline1[h];
    float *restrict sums0a = sums0 + h * n_samples;
    float *restrict sums1a = sums1 + h * n_samples;
    for (unsigned int i = 0; i < n_samples; i ++) {
      sums0a[i] = s0;
      sums1a[i] = s1;
    }
  }
}

static inline void tk_compressor_latent_sums_thread (
  uint64_t *restrict samples,
  unsigned int *restrict visibles,
  float *restrict alpha,
  float *restrict log_marg,
  float *restrict sums,
  unsigned int n_samples,
  unsigned int n_visible,
  unsigned int n_hidden,
  unsigned int bfirst,
  unsigned int blast
) {
  float *restrict sums0 = sums + 0 * n_hidden * n_samples;
  float *restrict sums1 = sums + 1 * n_hidden * n_samples;
  float *restrict lm00 = log_marg + 0 * n_hidden * n_visible;
  float *restrict lm01 = log_marg + 1 * n_hidden * n_visible;
  float *restrict lm10 = log_marg + 2 * n_hidden * n_visible;
  float *restrict lm11 = log_marg + 3 * n_hidden * n_visible;
  for (unsigned int b = bfirst; b <= blast; b ++) {
    uint64_t s = samples[b];
    unsigned int v = visibles[b];
    for (unsigned int h = 0; h < n_hidden; h ++) { // not vectorized, gather/scatter
      float *restrict sums0h = sums0 + h * n_samples;
      float *restrict sums1h = sums1 + h * n_samples;
      float *restrict lm00a = lm00 + h * n_visible;
      float *restrict lm01a = lm01 + h * n_visible;
      float *restrict lm10a = lm10 + h * n_visible;
      float *restrict lm11a = lm11 + h * n_visible;
      float *restrict aph = alpha + h * n_visible;
      sums0h[s] = sums0h[s] - aph[v] * lm00a[v] + aph[v] * lm10a[v]; // gather/scatter
      sums1h[s] = sums1h[s] - aph[v] * lm01a[v] + aph[v] * lm11a[v]; // gather/scatter
    }
  }
}

static inline void tk_compressor_latent_py_thread (
  float *restrict log_py,
  float *restrict log_pyx_unnorm,
  float *restrict sums,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  float *restrict sums0 = sums + 0 * n_hidden * n_samples;
  float *restrict sums1 = sums + 1 * n_hidden * n_samples;
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, not affine
    float lpy0v = log_py[h];
    float lpy1v = log1pf(-expf(lpy0v));
    float *restrict sums0h = sums0 + h * n_samples;
    float *restrict sums1h = sums1 + h * n_samples;
    float *restrict lpyx0 = log_pyx_unnorm + 0 * n_hidden * n_samples + h * n_samples;
    float *restrict lpyx1 = log_pyx_unnorm + 1 * n_hidden * n_samples + h * n_samples;
    for (unsigned int i = 0; i < n_samples; i ++) {
      lpyx0[i] = sums0h[i] + lpy0v;
      lpyx1[i] = sums1h[i] + lpy1v;
    }
  }
}

static inline void tk_compressor_latent_norm_thread (
  float *restrict log_z, // mis
  float *restrict pyx,
  float *restrict log_pyx_unnorm,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  float *restrict lpyx0 = log_pyx_unnorm + 0 * n_hidden * n_samples;
  float *restrict lpyx1 = log_pyx_unnorm + 1 * n_hidden * n_samples;
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++) {
    float a = lpyx0[i];
    float b = lpyx1[i];
    float max_ab  = (a > b) ? a : b;
    float sum_exp = expf(a - max_ab) + expf(b - max_ab);
    log_z[i] = max_ab + logf(sum_exp);
  }
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++)
    pyx[i] = expf(lpyx0[i] - log_z[i]);
}

static inline void tk_compressor_update_tc_thread (
  float *restrict log_z, // mis
  float *restrict tcs,
  unsigned int n_samples,
  unsigned int hfirst,
  unsigned int hlast
) {
  for (unsigned int h = hfirst; h <= hlast; h ++) { // not vectorized, unsupported outer form
    const float *restrict lz = log_z + h * n_samples;
    float sum = 0.0;
    for (unsigned int s = 0; s < n_samples; s ++)
      sum += lz[s];
    float tc = sum / (float) n_samples;
    tcs[h] = tc;
  }
}

static inline void tk_compressor_update_last_tc (
  float *restrict tcs,
  float *last_tc,
  unsigned int n_hidden
) {
  float sum0 = 0.0;
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
  float *restrict px,
  float *restrict entropy_x,
  unsigned int n_samples,
  unsigned int n_visible
) {
  for (unsigned int v = 0; v < n_visible; v ++)
    px[v] = 0;
  for (uint64_t c = 0; c < cardinality; c ++)
    px[visibles[c]] ++;
  for (unsigned int v = 0; v < n_visible; v ++)
    px[v] /= (float) n_samples;
  for (unsigned int v = 0; v < n_visible; v ++) {
    float entropy = 0;
    entropy -= px[v] * logf(px[v]);
    entropy -= (1 - px[v]) * logf(1 - px[v]);
    entropy_x[v] = entropy > 0 ? entropy : 1e-10;
  }
}

static inline void tk_compressor_init_alpha_thread (
  float *alpha,
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
  float *tcs,
  unsigned int hfirst,
  unsigned int hlast
) {
  for (unsigned int i = hfirst; i <= hlast; i ++)
    tcs[i] = 0.0;
}

static inline void tk_compressor_init_log_pyx_unnorm_thread (
  float *log_pyx_unnorm,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int hfirst,
  unsigned int hlast
) {
  float log_dim_hidden = -logf(2);
  float *restrict lpyx0 = log_pyx_unnorm + 0 * n_hidden * n_samples;
  float *restrict lpyx1 = log_pyx_unnorm + 1 * n_hidden * n_samples;
  for (unsigned int i = hfirst * n_samples; i < (hlast + 1) * n_samples; i ++) {
    lpyx0[i] = log_dim_hidden * (0.5 + fast_drand());
    lpyx1[i] = log_dim_hidden * (0.5 + fast_drand());
  }
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
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_MARGINALS:
        tk_compressor_marginals_thread(
          data->cardinality,
          data->C->samples,
          data->C->visibles,
          data->C->log_py,
          data->C->pyx,
          data->C->counts,
          data->C->log_marg,
          data->C->mis,
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
          data->C->alpha,
          data->C->baseline,
          data->C->log_marg,
          data->C->tcs,
          data->C->mis,
          data->C->maxmis,
          data->C->lam,
          data->C->n_visible,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_LATENT_BASELINE:
        tk_compressor_latent_baseline_thread(
          data->C->sums,
          data->C->baseline,
          data->n_samples,
          data->C->n_hidden,
          data->hfirst,
          data->hlast);
        break;
      case TK_CMP_LATENT_SUMS:
        tk_compressor_latent_sums_thread(
          data->C->samples,
          data->C->visibles,
          data->C->alpha,
          data->C->log_marg,
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
  C->mis = tk_realloc(L, C->mis, C->n_hidden * n_samples * sizeof(float));
  C->pyx = tk_realloc(L, C->pyx, C->n_hidden * n_samples * sizeof(float));
  C->log_pyx_unnorm = tk_realloc(L, C->log_pyx_unnorm, 2 * C->n_hidden * n_samples * sizeof(float));
  unsigned int len_sums = (2 * C->n_hidden * (n_samples > C->n_visible ? n_samples : C->n_visible));
  C->sums = tk_realloc(L, C->sums, len_sums * sizeof(float));
  tk_compressor_signal(
    TK_CMP_LATENT_BASELINE,
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
      float py0 = C->pyx[h * n_samples + s];
      if (py0 < 0.5)
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
  C->pyx = tk_malloc(L, C->n_hidden * n_samples * sizeof(float));
  C->log_pyx_unnorm = tk_malloc(L, 2 * C->n_hidden * n_samples * sizeof(float));
  C->sums = tk_malloc(L, 2 * C->n_hidden * n_samples * sizeof(float));
  unsigned int len_mis = C->n_hidden * C->n_visible;
  if (len_mis < (C->n_hidden * n_samples))
    len_mis = C->n_hidden * n_samples;
  C->mis = tk_malloc(L, len_mis * sizeof(float));
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
      TK_CMP_LATENT_BASELINE,
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
  C->tcs = tk_malloc(L, C->n_hidden * sizeof(float));
  C->alpha = tk_malloc(L, C->n_hidden * C->n_visible * sizeof(float));
  C->log_py = tk_malloc(L, C->n_hidden * sizeof(float));
  C->log_marg = tk_malloc(L, 2 * 2 * C->n_hidden * C->n_visible * sizeof(float));
  C->counts = tk_malloc(L, 2 * 2 * C->n_hidden * C->n_visible * sizeof(float));
  C->baseline = tk_malloc(L, 2 * C->n_hidden * sizeof(float));
  C->px = tk_malloc(L, C->n_visible * sizeof(float));
  C->entropy_x = tk_malloc(L, C->n_visible * sizeof(float));
  C->maxmis = tk_malloc(L, C->n_visible * sizeof(float));
  C->n_threads = n_threads;
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
  tk_lua_fwrite(L, &C->lam, sizeof(float), 1, fh);
  tk_lua_fwrite(L, &C->tmin, sizeof(float), 1, fh);
  tk_lua_fwrite(L, &C->ttc, sizeof(float), 1, fh);
  tk_lua_fwrite(L, C->alpha, sizeof(float), C->n_hidden * C->n_visible, fh);
  tk_lua_fwrite(L, C->log_py, sizeof(float), C->n_hidden, fh);
  tk_lua_fwrite(L, C->log_marg, sizeof(float), 2 * 2 * C->n_hidden * C->n_visible, fh);
  tk_lua_fwrite(L, C->baseline, sizeof(float), 2 * C->n_hidden, fh);
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
  tk_lua_fread(L, &C->lam, sizeof(float), 1, fh);
  tk_lua_fread(L, &C->tmin, sizeof(float), 1, fh);
  tk_lua_fread(L, &C->ttc, sizeof(float), 1, fh);
  C->alpha = tk_malloc(L, C->n_hidden * C->n_visible * sizeof(float));
  tk_lua_fread(L, C->alpha, sizeof(float), C->n_hidden * C->n_visible, fh);
  C->log_py = tk_malloc(L, C->n_hidden * sizeof(float));
  tk_lua_fread(L, C->log_py, sizeof(float), C->n_hidden, fh);
  C->log_marg = tk_malloc(L, 2 * 2 * C->n_hidden * C->n_visible * sizeof(float));
  tk_lua_fread(L, C->log_marg, sizeof(float), 2 * 2 * C->n_hidden * C->n_visible, fh);
  C->baseline = tk_malloc(L, 2 * C->n_hidden * sizeof(float));
  tk_lua_fread(L, C->baseline, sizeof(float), 2 * C->n_hidden, fh);
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
