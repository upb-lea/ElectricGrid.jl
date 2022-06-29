using Genie.Router
using NcConfigsController

route("/") do
  serve_static_file("welcome.html")
end

route("/nc", NcConfigsController.index)

route("/ncconfigs/search", NcConfigsController.search, named = :search_ncconfigs)