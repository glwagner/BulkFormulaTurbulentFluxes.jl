module SimilarityTheoryTurbulentFluxes

using Adapt
using Printf

using Thermodynamics: Liquid
using Thermodynamics: PhasePartition
using KernelAbstractions.Extras.LoopInfo: @unroll

import Thermodynamics as AtmosphericThermodynamics
import Thermodynamics.Parameters: molmass_ratio

include("seawater_saturation_specific_humidity.jl")
include("self_similar_boundary_layers.jl")
include("compute_turbulent_fluxes.jl")

end # module SimilarityTheoryTurbulentFluxes
