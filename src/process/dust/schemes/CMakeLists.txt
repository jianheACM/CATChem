# Dust Schemes CMakeLists.txt
# Generated on: 2025-08-05T10:07:05.835498

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - DustScheme_FENGSHA_Mod.F90
# - DustScheme_GINOUX_Mod.F90

# Scheme-specific configurations can be added here if needed

# FENGSHA scheme configuration
# Required meteorological fields for FENGSHA:
#   - {'name': 'IsLand', 'description': 'Land mask', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'is_land'}
#   - {'name': 'USTAR', 'description': 'Friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'friction_velocity'}
#   - {'name': 'LWI', 'description': 'Land-water index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'land_water_index'}
#   - {'name': 'GVF', 'description': 'Green vegetation fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'green_veg_fraction'}
#   - {'name': 'LAI', 'description': 'Leaf area index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'leaf_area_index'}
#   - {'name': 'FROCEAN', 'description': 'Ocean fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'ocean_fraction'}
#   - {'name': 'CLAYFRAC', 'description': 'Clay fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'clay_fraction'}
#   - {'name': 'SANDFRAC', 'description': 'Sand fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'sand_fraction'}
#   - {'name': 'FRSNO', 'description': 'Snow fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'snow_fraction'}
#   - {'name': 'RDRAG', 'description': 'Drag partition', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'drag_partition'}
#   - {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'}
#   - {'name': 'USTAR_THRESHOLD', 'description': 'Threshold friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'threshold_ustar'}

# GINOUX scheme configuration
# Required meteorological fields for GINOUX:
#   - {'name': 'FRLAKE', 'description': 'Lake fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'lake_fraction'}
#   - {'name': 'GWETTOP', 'description': 'Top soil wetness', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'top_soil_wetness'}
#   - {'name': 'U10M', 'description': '10m U wind', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'u_wind_10m'}
#   - {'name': 'V10M', 'description': '10m V wind', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'v_wind_10m'}
#   - {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'}


# Scheme validation targets
add_custom_target(validate_fengsha_scheme
    COMMAND ${CMAKE_COMMAND} -E echo "Validating fengsha scheme implementation"
    COMMENT "Validating fengsha scheme"
)
add_custom_target(validate_ginoux_scheme
    COMMAND ${CMAKE_COMMAND} -E echo "Validating ginoux scheme implementation"
    COMMENT "Validating ginoux scheme"
)

# Combined validation target
add_custom_target(validate_all_schemes
    DEPENDS
        validate_fengsha_scheme
        validate_ginoux_scheme
    COMMENT "Validating all dust schemes"
)