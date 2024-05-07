local env = {

  name = "santoku-bitmap",
  version = "0.0.8-1",
  variable_prefix = "TK_BITMAP",
  public = true,

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.204-1",
  },

  test = {
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
