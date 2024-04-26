struct ClasiusClapyeronSaturation end

#####
##### Models for the saturation specific humidity above pure water
#####
 
@inline function water_saturation_specific_humidity(::ClasiusClapyeronSaturation, ℂₐ, ρₛ, Tₛ)
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(ℂₐ, Tₛ, Liquid())
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(ℂₐ, Tₛ, ρₛ, p★)
    return q★
end

struct LargeYeagerSaturation{FT}
    c₁ :: FT
    c₂ :: FT
end

function LargeYeagerSaturation(FT=Float64; c₁ = 640380, c₂ = 5107.4)
    return LargeYeagerSaturation(convert(FT, c₁), convert(FT, c₂))
end

const LYS = LargeYeagerSaturation
@inline water_saturation_specific_humidity(lys::LYS, ℂₐ, ρₛ, Tₛ) = lys.c₁ * exp(-lys.c₂ / Tₛ) / ρₛ

#####
##### seawater saturation specific humidity
#####

@inline function seawater_saturation_specific_humidity(atmosphere_thermodynamics_parameters,
                                                       surface_temperature,
                                                       surface_salinity,
                                                       atmos_state,
                                                       water_mole_fraction,
                                                       water_vapor_saturation,
                                                       ::Liquid)

    ℂₐ = atmosphere_thermodynamics_parameters
    FT = eltype(ℂₐ)
    Tₛ = surface_temperature
    Sₛ = surface_salinity
    ρₛ = atmos_state.ρ # surface density -- should we extrapolate to obtain this?
    ρₛ = convert(FT, ρₛ)

    q★_H₂O = water_saturation_specific_humidity(water_vapor_saturation, ℂₐ, ρₛ, Tₛ)
    x_H₂O  = compute_water_mole_fraction(water_mole_fraction, Sₛ)

    # Return saturation specific humidity for salty seawater
    return q★_H₂O * x_H₂O
end

struct SalinityConstituent{FT}
    molar_mass :: FT
    mass_fraction :: FT
end

struct WaterMoleFraction{FT, C}
    water_molar_mass :: FT
    salinity_constituents :: C
end

function WaterMoleFraction(FT=Float64)
    water_molar_mass = convert(FT, 18.02)

    # TODO: find reference for these
    salinity_constituents = (
        chloride  = SalinityConstituent{FT}(35.45, 0.56),
        sodium    = SalinityConstituent{FT}(22.99, 0.31),
        sulfate   = SalinityConstituent{FT}(96.06, 0.08),
        magnesium = SalinityConstituent{FT}(24.31, 0.05),
    )

    return SeawaterComposition(water_molar_mass, salinity_constituents)
end

@inline compute_water_mole_fraction(x_H₂O::Number, S) = x_H₂O

@inline function compute_water_mole_fraction(wmf::WaterMoleFraction, S)
    # Molecular weights
    μ_H₂O = wmf.water_molar_mass

    # Salinity constituents: Cl, Na, SO₄, Mg
    μ_Cl  = wmf.salinity_constituents.chloride.molar_mass
    μ_Na  = wmf.salinity_constituents.sodium.molar_mass
    μ_SO₄ = wmf.salinity_constituents.sulfate.molar_mass
    μ_Mg  = wmf.salinity_constituents.magnesium.molar_mass

    # Salinity constituent fractions
    ϵ_Cl  = wmf.salinity_constituents.chloride.mass_fraction
    ϵ_Na  = wmf.salinity_constituents.sodium.mass_fraction
    ϵ_SO₄ = wmf.salinity_constituents.sulfate.mass_fraction
    ϵ_Mg  = wmf.salinity_constituents.magnesium.mass_fraction

    α = μ_H₂O * (ϵ_Cl/μ_Cl + ϵ_Na/μ_Na  + ϵ_SO₄/μ_SO₄ + ϵ_Mg/μ_Mg)

    # TODO: this assumes that `S` has units of g / kg.
    # Should we express the concept of "ocean_salinity_units"?
    ς = S / 1000 # convert g/kg to concentration

    return (1 - ς) / (1 - ς + α * ς)
end

