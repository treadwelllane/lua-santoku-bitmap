local env = {

  name = "santoku-bitmap",
  version = "0.0.22-1",
  variable_prefix = "TK_BITMAP",
  public = true,

  -- -Wsign-conversion -Wsign-compare would be nice, but need to fix all the warning spots
  cflags = { "-march=native", "-O3", "-ffast-math", "-Wall", "-Wextra", "-Wno-sign-compare", "-Wstrict-overflow", "-Wpointer-sign", "-Wno-unused-parameter" },
  ldflags = { "-march=native", "-O3", "-lm", },

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.245-1",
  },

  test = {
    cflags = { "-fopt-info-vec=optimize.txt", "-fopt-info-vec-missed=optimize.txt", "-g3" },
    ldflags = { "-fopt-info-vec=optimize.txt", "-fopt-info-vec-missed=optimize.txt", "-g3" },
    dependencies = {
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
