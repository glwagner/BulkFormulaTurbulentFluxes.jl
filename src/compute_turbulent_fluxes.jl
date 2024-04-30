@inline compute_similarity_theory_fluxes(turbulent_fluxes, atmos_state, surface_state) =
    compute_similarity_theory_fluxes(turbulent_fluxes.roughness_lengths, turbulent_fluxes, atmos_state, surface_state)

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct BoundaryLayerScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

#####
##### Fixed-point iteration for roughness length
#####

struct SimilarityFunction{FT, C}
    a :: FT
    b :: FT
    c :: C
end

@inline function (ψ::SimilarityFunction)(Ri)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    Ri⁻ = min(zero(Ri), Ri)
    ϕ⁻¹ = (1 - b * Ri⁻)^c
    ψ_unstable = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - (4 * atan(ϕ⁻¹) + π) / 2

    ψ_stable = - a * Ri

    return ifelse(Ri < 0, ψ_unstable, ψ_stable)
end

struct OneQuarter end
struct OneHalf end

import Base: ^
@inline ^(x, ::OneQuarter) = sqrt(sqrt(x))
@inline ^(x, ::OneHalf) = sqrt(x)

@inline function bulk_factor(ψ, h, ℓ, Ri)
    L★ = h / Ri
    χ⁻¹ = log(h / ℓ) - ψ(Ri) + ψ(ℓ / L★)
    return 1 / χ⁻¹
end

@inline function buoyancy_scale(θ★, q★, 𝒬, ℂ, g)
    𝒯₀ = AtmosphericThermodynamics.virtual_temperature(ℂ, 𝒬)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬)

    ε = AtmosphericThermodynamics.Parameters.molmass_ratio(ℂ)
    δ = ε - 1

    b★ = g / 𝒯₀ * (θ★ * (1 + δ * q₀) + δ * θ₀ * q★)

    return b★
end

default_initial_similarity_scales() = (momentum=1e-3, temperature=1e-3, water_vapor=1e-3)

@inline function compute_similarity_theory_fluxes(roughness_lengths,
                                                  surface_state,
                                                  atmos_state,
                                                  thermodynamics_parameters,
                                                  gravitational_acceleration,
                                                  von_karman_constant,
                                                  initial_similarity_scales = default_initial_similarity_scales)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, atmos_state, surface_state)
    differences = (; u=Δu, v=Δv, θ=Δθ, q=Δq, h=Δh)

    # Solve for the characteristic scales u★, θ★, q★, and thus for fluxes.
    Σ★ = initial_similarity_scales
 
    @unroll for iter = 1:10
        Σ★ = refine_characteristic_scales(Σ★,
                                          roughness_lengths, 
                                          surface_state,
                                          differences,
                                          thermodynamics_parameters,
                                          gravitational_acceleration,
                                          von_karman_constant)
    end

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    # u★² ≡ sqrt(τx² + τy²)
    τx = - u★^2 * Δu / sqrt(Δu^2 + Δv^2)
    τy = - u★^2 * Δv / sqrt(Δu^2 + Δv^2)

    𝒬ₐ = atmos_state.ts
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    fluxes = (;
        water_vapor   = - ρₐ * u★ * q★,
        sensible_heat = - ρₐ * cₚ * u★ * θ★,
        latent_heat   = - ρₐ * u★ * q★ * ℰv,
        x_momentum    = + ρₐ * τx,
        y_momentum    = + ρₐ * τy,
    )

    return fluxes
end

@inline compute_roughness_length(ℓ::Number, Σ★) = ℓ
@inline compute_roughness_length(ℓ, Σ★) = ℓ(Σ★)

@inline function refine_characteristic_scales(estimated_characteristic_scales,
                                              roughness_lengths,
                                              surface_state,
                                              differences,
                                              thermodynamics_parameters,
                                              gravitational_acceleration,
                                              von_karman_constant)

    # "initial" scales because we will recompute them
    u★ = estimated_characteristic_scales.momentum
    θ★ = estimated_characteristic_scales.temperature
    q★ = estimated_characteristic_scales.water_vapor
    Σ★ = estimated_characteristic_scales

    # Extract roughness lengths
    ℓu = roughness_lengths.momentum
    ℓθ = roughness_lengths.temperature
    ℓq = roughness_lengths.water_vapor

    ℓu₀ = compute_roughness_length(ℓu, Σ★)
    ℓθ₀ = compute_roughness_length(ℓθ, Σ★)
    ℓq₀ = compute_roughness_length(ℓq, Σ★)

    # Compute flux Richardson number
    h = differences.h
    ϰ = von_karman_constant

    ℂ = thermodynamics_parameters
    g = gravitational_acceleration
    𝒬ₒ = surface_state.ts # thermodyanmic state
    b★ = buoyancy_scale(θ★, q★, 𝒬ₒ, ℂ, g)
    Riₕ = - ϰ * h * b★ / u★^2
    Riₕ = ifelse(isnan(Riₕ), zero(Riₕ), Riₕ) 

    # Compute similarity functions
    ψu = SimilarityFunction(4.7, 15.0, OneQuarter())
    ψc = SimilarityFunction(6.35, 9.0, OneHalf())

    χu = bulk_factor(ψu, h, ℓu₀, Riₕ)
    χθ = bulk_factor(ψc, h, ℓθ₀, Riₕ)
    χq = bulk_factor(ψc, h, ℓq₀, Riₕ)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    u★ = ϰ * χu * sqrt(Δu^2 + Δv^2)
    θ★ = ϰ * χθ * Δθ
    q★ = ϰ * χq * Δq

    return (momentum=u★, temperature=θ★, water_vapor=q★)
end

struct GravityWaveRoughnessLength{FT}
    gravitational_acceleration :: FT
    air_kinematic_viscosity :: FT
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
end

function GravityWaveRoughnessLength(FT=Float64;
                                    gravitational_acceleration = default_gravitational_acceleration,
                                    air_kinematic_viscosity = 1.5e-5,
                                    gravity_wave_parameter = 0.011,
                                    laminar_parameter = 0.11)

    return GravityWaveRoughnessLength(convert(FT, gravitational_acceleration),
                                      convert(FT, air_kinematic_viscosity),
                                      convert(FT, gravity_wave_parameter),
                                      convert(FT, laminar_parameter))
end

@inline function compute_roughness_length(ℓ::GravityWaveRoughnessLength, Σ★)
    u★ = Σ★.momentum
    g = ℓ.gravitational_acceleration
    ν = ℓ.air_kinematic_viscosity
    α = ℓ.gravity_wave_parameter
    β = ℓ.laminar_parameter

    return α * u★^2 / g + β * ν / u★
end

function default_roughness_lengths(FT=Float64)
    momentum    = GravityWaveRoughnessLength(FT)
    temperature = convert(FT, 1e-4)
    water_vapor = convert(FT, 1e-4)
    return BoundaryLayerScales(momentum, temperature, water_vapor)
end

