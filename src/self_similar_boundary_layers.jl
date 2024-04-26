#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SelfSimilarBoundaryLayer{FT, ΔU, UF, TP, S, W, R, F}
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    bulk_velocity_scale :: ΔU
    similarity_functions :: UF
    thermodynamics_parameters :: TP
    water_vapor_saturation :: S
    water_mole_fraction :: W
    roughness_lengths :: R
end

@inline molmass_ratio(fluxes::STTF) = molmass_ratio(fluxes.thermodynamics_parameters)

@inline universal_func_type(fluxes::STTF{<:Any, <:Any, <:BusingerParams}) = BusingerType()

Adapt.adapt_structure(to, fluxes::STTF) = SelfSimilarBoundaryLayer(adapt(to, fluxes.gravitational_acceleration),
                                                                          adapt(to, fluxes.von_karman_constant),
                                                                          nothing, # adapt(to, fluxes.bulk_velocity_scale),
                                                                          adapt(to, fluxes.similarity_functions),
                                                                          adapt(to, fluxes.thermodynamics_parameters),
                                                                          nothing, #adapt(to, fluxes.water_vapor_saturation),
                                                                          nothing, #adapt(to, fluxes.water_mole_fraction),
                                                                          adapt(to, fluxes.roughness_lengths))

Base.summary(::SelfSimilarBoundaryLayer{FT}) where FT = "SelfSimilarBoundaryLayer{$FT}"

function Base.show(io::IO, fluxes::SelfSimilarBoundaryLayer)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ",   fluxes.gravitational_acceleration, '\n',
          "├── von_karman_constant: ",          fluxes.von_karman_constant, '\n',
          "├── bulk_velocity_scale: ",          summary(fluxes.bulk_velocity_scale), '\n',
          "├── similarity_function: ",          summary(fluxes.similarity_function), '\n',
          "├── water_mole_fraction: ",          summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",       summary(fluxes.water_vapor_saturation), '\n',
          "└── thermodynamics_parameters: ",    summary(fluxes.thermodynamics_parameters))
end

# const PATP = PrescribedAtmosphereThermodynamicsParameters

function SelfSimilarBoundaryLayer(FT::DataType = Float64;
                                  gravitational_acceleration = default_gravitational_acceleration,
                                  bulk_velocity_scale = nothing,
                                  von_karman_constant = convert(FT, 0.4),
                                  similarity_functions = businger_similarity_functions(FT),
                                  thermodynamics_parameters = nothing,
                                  water_vapor_saturation = ClasiusClapyeronSaturation(),
                                  water_mole_fraction = convert(FT, 0.98),
                                  roughness_lengths = default_roughness_lengths(FT))

    return SelfSimilarBoundaryLayer(convert(FT, gravitational_acceleration),
                                           convert(FT, von_karman_constant),
                                           bulk_velocity_scale,
                                           similarity_functions,
                                           thermodynamics_parameters,
                                           water_vapor_saturation,
                                           water_mole_fraction,
                                           roughness_lengths)
end

