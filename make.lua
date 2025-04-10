local env = {

  name = "santoku-bitmap",
  version = "0.0.50-1",
  variable_prefix = "TK_BITMAP",
  public = true,

  cflags = { "-march=native", "-O3", "-Wall", "-Wextra", "-Wno-sign-conversion", "-Wsign-compare", "-Wstrict-overflow", "-Wpointer-sign", },
  ldflags = { "-march=native", "-O3", "-lm", "-lpthread", "-lnuma" },

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.247-1",
  },

  test = {
    cflags = { "-g3" },
    ldflags = { "-g3" },
    dependencies = {
      "santoku-fs == 0.0.34-1",
      "santoku-matrix == 0.0.13-1",
      "luacov == 0.15.0-1",
    }
  },

}

env.homepage = "https://github.com/treadwelllane/lua-" .. env.name
env.tarball = env.name .. "-" .. env.version .. ".tar.gz"
env.download = env.homepage .. "/releases/download/" .. env.version .. "/" .. env.tarball

return {
  type = "lib",
  env = env,
}
