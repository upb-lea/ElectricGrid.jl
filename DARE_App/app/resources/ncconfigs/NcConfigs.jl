module NcConfigs

import SearchLight: AbstractModel, DbId
import Base: @kwdef

export NcConfig

@kwdef mutable struct NcConfig <: AbstractModel
  id::DbId = DbId()
  name::String = ""
  cm::String = ""
  parameters::String = ""
end

end
