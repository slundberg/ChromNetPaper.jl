module ChromNetPaper

globalDataDir = normpath(joinpath(dirname(Base.source_path()), "..", "data"))

include("evaluate.jl")
include("groupgm.jl")
include("drive.jl")

end # module
