
struct ConstantExchangeCoefficients{FT, ΔU, TP, S, W}
    momentum_exchange_coefficient :: FT
    heat_exchange_coefficient :: FT
    water_vapor_exchange_coefficient :: FT
    bulk_velocity_scale :: ΔU
    thermodynamics_parameters :: TP
    water_vapor_saturation :: S
    water_mole_fraction :: W
end

function ConstantExchangeCoefficients(FT::DataType = Float64;
                                      momentum_exchange_coefficient = 1e-3,
                                      heat_exchange_coefficient = 1e-3,
                                      water_vapor_exchange_coefficient = 1e-3,
                                      bulk_velocity_scale = nothing,
                                      thermodynamics_parameters = nothing,
                                      water_vapor_saturation = ClasiusClapyeronSaturation(),
                                      water_mole_fraction = convert(FT, 0.98))

    return ConstantExchangeCoefficients(convert(FT, momentum_exchange_coefficient), 
                                        convert(FT, heat_exchange_coefficient),
                                        convert(FT, water_vapor_exchange_coefficient),
                                        bulk_velocity_scale,
                                        thermodynamics_parameters,
                                        water_vapor_saturation,
                                        water_mole_fraction = convert(FT, 0.98))
end

function compute_turbulent_fluxes(constant_coefficients::ConstantExchangeCoeffients,
                                  surface_state,
                                  atmos_state,
                                  thermodynamics_parameters)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, atmos_state, surface_state)

    Cm = constant_coefficients.momentum_exchange_coefficient
    Ch = constant_coefficients.heat_exchange_coefficient
    Cq = constant_coefficients.water_vapor_exchange_coefficient

    𝒬ₐ = atmos_state.ts
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    τx = - Cm * Δu * sqrt(Δu^2 + Δv^2)
    τy = - Cm * Δv * sqrt(Δu^2 + Δv^2)
    Jᶿ = - Ch * Δθ
    F  = - Cq * Δq

    fluxes = (;
        water_vapor   = - ρₐ * F,
        sensible_heat = - ρₐ * cₚ * Jᶿ,
        latent_heat   = - ρₐ * ℰv * F,
        x_momentum    = + ρₐ * τx,
        y_momentum    = + ρₐ * τy,
   )

    return fluxes
end

