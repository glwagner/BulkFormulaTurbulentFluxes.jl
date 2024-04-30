module SimilarityTheoryTurbulentFluxes

export compute_turbulent_fluxes

using Adapt
using Printf

using Thermodynamics: Liquid
using Thermodynamics: PhasePartition
using KernelAbstractions.Extras.LoopInfo: @unroll

import Thermodynamics as AtmosphericThermodynamics
import Thermodynamics.Parameters: molmass_ratio

@inline function state_differences(ℂ, 𝒰₁, 𝒰₀)
    z₁ = 𝒰₁.z
    z₀ = 𝒰₀.z
    Δh = z₁ - z₀

    U₁ = 𝒰₁.u
    U₀ = 𝒰₀.u

    @inbounds begin
        Δu = U₁[1] - U₀[1]
        Δv = U₁[2] - U₀[2]
    end

    # Thermodynamic state
    𝒬₁ = 𝒰₁.ts
    𝒬₀ = 𝒰₀.ts

    θ₁ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₁)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)
    Δθ = θ₁ - θ₀

    q₁ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₁)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₀)
    Δq = q₁ - q₀

    return Δh, Δu, Δv, Δθ, Δq
end

include("seawater_saturation_specific_humidity.jl")
include("constant_exchange_coefficients.jl")

# include("self_similar_boundary_layers.jl")
# include("compute_turbulent_fluxes.jl")

end # module SimilarityTheoryTurbulentFluxes
