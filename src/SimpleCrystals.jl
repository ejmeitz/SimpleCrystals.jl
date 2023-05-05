module SimpleCrystals

using Reexport
using PeriodicTable

@reexport using AtomsBase
@reexport using Unitful
@reexport using StaticArrays

include("atom.jl")
include("BravaisLattice.jl")
include("Crystal.jl")
include("FileWriter.jl")

include("bravais/2D_bravais.jl")
include("bravais/3D_bravais.jl")
include("non_bravais/2D_other.jl")
include("non_bravais/3D_other.jl")

end