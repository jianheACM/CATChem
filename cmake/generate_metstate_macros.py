#!/usr/bin/env python3
"""
Generate Fortran macros for MetStateType fields.

This script parses a Fortran source file containing the MetStateType definition and generates Fortran include files
for field accessors, allocation, and deallocation macros. The generated macros are used to implement field-by-field
allocation and access in a modular, robust way.

Usage
-----
Run from the command line:
    generate_metstate_macros.py accessor <input_f90> <output_inc>
    generate_metstate_macros.py allocate <input_f90> <output_inc>
    generate_metstate_macros.py deallocate <input_f90> <output_inc>
    generate_metstate_macros.py column_accessor <input_f90> <output_inc>
    generate_metstate_macros.py 2d_scalar_accessor <input_f90> <output_inc>
    generate_metstate_macros.py scalar_accessor <input_f90> <output_inc>
    generate_metstate_macros.py virtualcolumn_populate <input_f90> <output_inc>
    generate_metstate_macros.py virtualmet_populate <input_f90> <output_inc>

Parameters
----------
mode : str
    Which macro to generate. One of: accessor, allocate, deallocate, column_accessor, 2d_scalar_accessor, scalar_accessor, virtualcolumn_populate, virtualmet_populate.
input_f90 : str
    Path to the Fortran source file containing the MetStateType definition.
output_inc : str
    Path to the output .inc file to write.

Examples
--------
$ python generate_metstate_macros.py allocate ../src/core/metstate_mod.F90 metstate_allocate_fields.inc

"""
import re
import sys

def parse_metstate_type(filename):
    """
    Parse the MetStateType definition in a Fortran source file.

    Parameters
    ----------
    filename : str
        Path to the Fortran source file containing MetStateType.

    Returns
    -------
    fields : list of tuple
        List of (name, rank, dims) for each allocatable REAL(fp) field.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    in_type = False
    fields = []  # List of (name, rank, dims)
    for line in lines:
        if 'TYPE, PUBLIC :: MetStateType' in line:
            in_type = True
        elif in_type and 'end type' in line.lower():
            break
        elif in_type:
            # Match real(fp), allocatable :: name(dimensions)
            m = re.match(r'\s*REAL\(fp\),\s*ALLOCATABLE\s*::\s*(\w+)\s*\(([^)]*)\)', line, re.IGNORECASE)
            if m:
                name = m.group(1)
                dims = m.group(2)
                rank = dims.count(',') + 1
                fields.append((name, rank, dims))

    return fields

def write_accessor(fields, output_file):
    """
    Write a macro for field accessors (all 2D/3D fields).

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    # Deduplicate by case-insensitive field name
    unique_fields = {}
    for name, rank, dims in fields:
        key = name.lower()
        if key not in unique_fields:
            unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case(" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) then\n")
            if rank == 3:
                f.write(f"      column_ptr => this%{name}(col_i, col_j, :)\n")
            elif rank == 2:
                f.write(f"      column_ptr => null()\n")
            f.write(f"   endif\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_allocate(fields, output_file):
    """
    Write a macro for allocation statements for MetStateType fields.
    Skips scalar (rank-0) fields. Adds a case for 'ALL'.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    with open(output_file, 'w') as f:
        # ALL case
        f.write("case ('ALL', 'all')\n")
        for name, rank, dims in fields:
            if rank == 0:
                continue
            if name.lower() in ("soilm"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoil))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frsoil"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoiltype))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frlanduse", "frlai", "frz0"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nSURFTYPE))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            else:
                if rank == 2:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz))\n")
                else:
                    raise ValueError(f"Unsupported rank {rank} for field {name}")
        # Individual field cases
        for name, rank, dims in fields:
            if rank == 0:
                continue  # Skip scalars
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            if name.lower() in ("soilm"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoil))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frsoil"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoiltype))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frlanduse", "frlai", "frz0"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nSURFTYPE))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            else:
                if rank == 2:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz))\n")
                else:
                    raise ValueError(f"Unsupported rank {rank} for field {name}")
        f.write("case default\n")
        f.write("  ! No allocation for unknown field\n")

def write_deallocate(fields, output_file):
    """
    Write a macro for deallocation statements for MetStateType fields.
    Skips scalar (rank-0) fields. Adds a case for 'ALL'.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    with open(output_file, 'w') as f:
        # ALL case
        f.write("case ('ALL', 'all')\n")
        for name, rank, dims in fields:
            if rank == 0:
                continue
            f.write(f"  if (allocated(this%{name})) deallocate(this%{name})\n")
        # Individual field cases
        for name, rank, dims in fields:
            if rank == 0:
                continue  # Skip scalars
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"  if (allocated(this%{name})) deallocate(this%{name})\n")
        f.write("case default\n")
        f.write("  ! No deallocation for unknown field\n")

def write_column_accessor(fields, output_file):
    """
    Write a macro for 3D column field accessors.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    # Only 3D fields
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 3:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) column_ptr => this%{name}(col_i, col_j, :)\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_2d_scalar_accessor(fields, output_file):
    """
    Write a macro for 2D scalar field accessors.

    Parameters
    ----------
    fields : list of (name, rank, dims)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    # Only 2D fields
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 2:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) scalar_val = this%{name}(col_i, col_j)\n")
        f.write("case default\n")
        f.write("   scalar_val = 0.0_fp\n")

def write_virtualmet_populate(fields, output_file):
    """
    Write a macro for populating VirtualMetType with MetState field pointers.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    
    # Define the 3D atmospheric fields (vertical levels)
    wanted_3d_fields = [
        'T', 'U', 'V', 'QV', 'PMID', 'AIRDEN', 'THETA', 'RH', 'OMEGA', 'CLDF',
        'Z', 'ZMID', 'BXHEIGHT', 'DELP', 'AIRVOL', 'MAIRDEN', 'PEDGE', 
        'PMID_DRY', 'DELP_DRY', 'PEDGE_DRY', 'QI', 'QL', 'CMFMC', 'DTRAIN',
        'F_OF_PBL', 'DQRCU', 'DQRLSAN', 'PFICU', 'PFILSAN', 'PFLCU', 
        'PFLLSAN', 'TAUCLI', 'TAUCLW', 'AIRNUMDEN', 'AVGW', 'SPHU',
        'TV', 'DAIRMASS', 'F_UNDER_PBLTOP'
    ]
    
    # Define categorical 3D fields (different third dimension)
    categorical_3d_fields = [
        'SOILM',       # (nx,ny,nsoil) - soil layers
        'FRLANDUSE',   # (nx,ny,nlanduse) - land use categories
        'FRSOIL',      # (nx,ny,nsoil) - soil type fractions
        'FRLAI',       # (nx,ny,nlanduse) - LAI per land use category
        'FRZ0'         # (nx,ny,nlanduse) - roughness per land use category
    ]
    
    # Define the 2D fields we want to include
    wanted_2d_fields = [
        'T2M', 'PS', 'USTAR', 'U10M', 'V10M', 'PBLH', 'Z0', 'LAI', 'FRLAND',
        'FROCEAN', 'FRSEAICE', 'AREA_M2', 'SST', 'TSKIN', 'TS', 'QV2M', 'SLP',
        'PHIS', 'SUNCOS', 'SWGDN', 'HFLUX', 'EFLUX', 'PRECCON', 'PRECLSC',
        'PRECANV', 'TO3', 'TROPP', 'TropHt'
    ]
    
    with open(output_file, 'w') as f:
        f.write("! Generated macro for populating VirtualMetType with MetState field pointers\n")
        f.write("! This macro should be included in the populate_virtual_column subroutine\n\n")
        
        # Populate 3D atmospheric field pointers
        f.write("! Populate 3D atmospheric field pointers (vertical levels: nlev)\n")
        for name, rank, dims in fields:
            if rank == 3 and name in wanted_3d_fields:
                f.write(f"call this%met_state%get_field_ptr('{name}', grid_i, grid_j, &\n")
                f.write(f"                                 col_ptr=column_ptr, rc=field_rc)\n")
                f.write(f"if (field_rc == 0 .and. associated(column_ptr)) then\n")
                f.write(f"   virtual_col%met%{name} => column_ptr\n")
                f.write(f"end if\n\n")
        
        # Populate categorical 3D field pointers  
        f.write("! Populate 3D categorical field pointers (NOT vertical levels!)\n")
        f.write("! WARNING: These have different dimensions - do not iterate with nlev!\n")
        for name, rank, dims in fields:
            if rank == 3 and name in categorical_3d_fields:
                f.write(f"call this%met_state%get_field_ptr('{name}', grid_i, grid_j, &\n")
                f.write(f"                                 col_ptr=column_ptr, rc=field_rc)\n")
                f.write(f"if (field_rc == 0 .and. associated(column_ptr)) then\n")
                f.write(f"   virtual_col%met%{name} => column_ptr\n")
                f.write(f"end if\n\n")
        
        # Populate 2D scalar fields
        f.write("! Populate 2D scalar fields\n")
        for name, rank, dims in fields:
            if rank == 2 and name in wanted_2d_fields:
                f.write(f"call this%met_state%get_field_ptr('{name}', grid_i, grid_j, &\n")
                f.write(f"                                 scalar_val=scalar_val, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val\n")
                f.write(f"end if\n\n")
        
        f.write("! Note: Categorical 3D fields (SOILM, FRLANDUSE, FRSOIL, FRLAI, FRZ0)\n")
        f.write("! are excluded as they require specific category indices.\n")
        f.write("! Scalar fields are not included as they are global, not column-specific.\n")

def write_virtualmet_type(fields, output_file):
    """
    Write the VirtualMetType definition macro based on MetState field definitions.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    
    # Define the 3D atmospheric fields (vertical levels)
    wanted_3d_fields = [
        'T', 'U', 'V', 'QV', 'PMID', 'AIRDEN', 'THETA', 'RH', 'OMEGA', 'CLDF',
        'Z', 'ZMID', 'BXHEIGHT', 'DELP', 'AIRVOL', 'MAIRDEN', 'PEDGE', 
        'PMID_DRY', 'DELP_DRY', 'PEDGE_DRY', 'QI', 'QL', 'CMFMC', 'DTRAIN',
        'F_OF_PBL', 'DQRCU', 'DQRLSAN', 'PFICU', 'PFILSAN', 'PFLCU', 
        'PFLLSAN', 'TAUCLI', 'TAUCLW', 'AIRNUMDEN', 'AVGW', 'SPHU',
        'TV', 'DAIRMASS', 'F_UNDER_PBLTOP'
    ]
    
    # Define categorical 3D fields (different third dimension - NOT vertical levels)
    categorical_3d_fields = [
        'SOILM',       # (nx,ny,nsoil) - soil layers
        'FRLANDUSE',   # (nx,ny,nlanduse) - land use categories
        'FRSOIL',      # (nx,ny,nsoil) - soil type fractions
        'FRLAI',       # (nx,ny,nlanduse) - LAI per land use category
        'FRZ0'         # (nx,ny,nlanduse) - roughness per land use category
    ]
    
    # Define the 2D fields we want to include
    wanted_2d_fields = [
        'T2M', 'PS', 'USTAR', 'U10M', 'V10M', 'PBLH', 'Z0', 'LAI', 'FRLAND',
        'FROCEAN', 'FRSEAICE', 'AREA_M2', 'SST', 'TSKIN', 'TS', 'QV2M', 'SLP',
        'PHIS', 'SUNCOS', 'SWGDN', 'HFLUX', 'EFLUX', 'PRECCON', 'PRECLSC',
        'PRECANV', 'TO3', 'TROPP', 'TropHt'
    ]
    
    # Field descriptions for better documentation
    field_descriptions = {
        'T': 'Temperature [K]',
        'U': 'U-wind [m/s]',
        'V': 'V-wind [m/s]',
        'QV': 'Water vapor [kg/kg]',
        'PMID': 'Mid-level pressure [Pa]',
        'AIRDEN': 'Air density [kg/m3]',
        'THETA': 'Potential temperature [K]',
        'RH': 'Relative humidity [%]',
        'OMEGA': 'Vertical velocity [Pa/s]',
        'CLDF': 'Cloud fraction [1]',
        'Z': 'Geopotential height [m]',
        'ZMID': 'Mid-level height [m]',
        'BXHEIGHT': 'Box height [m]',
        'DELP': 'Pressure thickness [Pa]',
        'AIRVOL': 'Air volume [m3]',
        'MAIRDEN': 'Moist air density [kg/m3]',
        'PEDGE': 'Edge pressure [Pa]',
        'PMID_DRY': 'Dry mid pressure [Pa]',
        'DELP_DRY': 'Dry pressure thickness [Pa]',
        'PEDGE_DRY': 'Dry edge pressure [Pa]',
        'QI': 'Ice water [kg/kg]',
        'QL': 'Liquid water [kg/kg]',
        'CMFMC': 'Cloud mass flux [kg/m2/s]',
        'DTRAIN': 'Detrainment flux [kg/m2/s]',
        'F_OF_PBL': 'Fraction in PBL [1]',
        'DQRCU': 'Conv precip rate [kg/kg/s]',
        'DQRLSAN': 'LS precip rate [kg/kg/s]',
        'PFICU': 'Conv ice precip flux [kg/m2/s]',
        'PFILSAN': 'LS ice precip flux [kg/m2/s]',
        'PFLCU': 'Conv liq precip flux [kg/m2/s]',
        'PFLLSAN': 'LS liq precip flux [kg/m2/s]',
        'TAUCLI': 'Ice cloud optical depth [1]',
        'TAUCLW': 'Water cloud optical depth [1]',
        'AIRNUMDEN': 'Air number density [molec/cm3]',
        'AVGW': 'Average molecular weight [g/mol]',
        'SPHU': 'Specific humidity [kg/kg]',
        'TV': 'Virtual temperature [K]',
        'DAIRMASS': 'Dry air mass [kg]',
        'F_UNDER_PBLTOP': 'Fraction under PBL top [1]',
        'SOILM': 'Volumetric soil moisture [m3/m3] (nsoil layers)',
        'FRLANDUSE': 'Fractional land use [1] (nlanduse categories)', 
        'FRSOIL': 'Fractional soil [1] (nsoil categories)',
        'FRLAI': 'LAI per land use type [m2/m2] (nlanduse categories)',
        'FRZ0': 'Roughness per land use [m] (nlanduse categories)',
        'T2M': '2m temperature [K]',
        'PS': 'Surface pressure [Pa]',
        'USTAR': 'Friction velocity [m/s]',
        'U10M': '10m U-wind [m/s]',
        'V10M': '10m V-wind [m/s]',
        'PBLH': 'PBL height [m]',
        'Z0': 'Roughness length [m]',
        'LAI': 'Leaf area index [m2/m2]',
        'FRLAND': 'Land fraction [1]',
        'FROCEAN': 'Ocean fraction [1]',
        'FRSEAICE': 'Sea ice fraction [1]',
        'AREA_M2': 'Grid area [m2]',
        'SST': 'Sea surface temperature [K]',
        'TSKIN': 'Skin temperature [K]',
        'TS': 'Surface temperature [K]',
        'QV2M': '2m water vapor [kg/kg]',
        'SLP': 'Sea level pressure [Pa]',
        'PHIS': 'Surface geopotential [m2/s2]',
        'SUNCOS': 'Cosine solar zenith angle [1]',
        'SWGDN': 'Downward SW radiation [W/m2]',
        'HFLUX': 'Sensible heat flux [W/m2]',
        'EFLUX': 'Latent heat flux [W/m2]',
        'PRECCON': 'Convective precipitation [kg/m2/s]',
        'PRECLSC': 'Large-scale precipitation [kg/m2/s]',
        'PRECANV': 'Anvil precipitation [kg/m2/s]',
        'TO3': 'Total ozone [DU]',
        'TROPP': 'Tropopause pressure [Pa]',
        'TropHt': 'Tropopause height [m]'
    }
    
    with open(output_file, 'w') as f:
        f.write("! Generated VirtualMetType definition based on MetState field definitions\n")
        f.write("! This macro should be included in the VirtualColumn_Mod.F90 type definition\n\n")
        
        f.write("   !> \\brief Virtual meteorological data container with direct pointers\n")
        f.write("   !! \\details Contains pointers to meteorological fields for a single column.\n")
        f.write("   !! Pointers are set directly from MetState accessor functions, eliminating\n")
        f.write("   !! data copying and providing efficient field access.\n")
        f.write("   type :: VirtualMetType\n")
        
        # Write 3D atmospheric field pointers (vertical levels)
        f.write("      ! 3D atmospheric fields (vertical profiles: nlev) - pointers to MetState data\n")
        for name, rank, dims in fields:
            if rank == 3 and name in wanted_3d_fields:
                comment = field_descriptions.get(name, f"{name} field")
                f.write(f"      real(fp), pointer :: {name}(:) => null()  !< {comment}\n")
        
        # Write categorical 3D field pointers (different dimensions)
        f.write("\n      ! 3D categorical fields (NOT vertical levels) - pointers to MetState data\n")
        f.write("      ! WARNING: These have different third dimensions (nsoil, nlanduse, etc.)\n")
        f.write("      ! Do NOT iterate with nlev - use appropriate category counts!\n")
        for name, rank, dims in fields:
            if rank == 3 and name in categorical_3d_fields:
                comment = field_descriptions.get(name, f"{name} field")
                f.write(f"      real(fp), pointer :: {name}(:) => null()  !< {comment}\n")
        
        f.write("\n      ! 2D surface fields (scalars for the column) - values copied from MetState\n")
        for name, rank, dims in fields:
            if rank == 2 and name in wanted_2d_fields:
                comment = field_descriptions.get(name, f"{name} field")
                f.write(f"      real(fp) :: {name} = 0.0_fp  !< {comment}\n")
        
        f.write("\n   contains\n")
        f.write("      procedure :: cleanup => virtual_met_cleanup\n")
        f.write("   end type VirtualMetType\n\n")

def write_virtualmet_cleanup(fields, output_file):
    """
    Write the VirtualMetType cleanup procedure macro.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    
    # Define the 3D atmospheric fields (vertical levels)
    wanted_3d_fields = [
        'T', 'U', 'V', 'QV', 'PMID', 'AIRDEN', 'THETA', 'RH', 'OMEGA', 'CLDF',
        'Z', 'ZMID', 'BXHEIGHT', 'DELP', 'AIRVOL', 'MAIRDEN', 'PEDGE', 
        'PMID_DRY', 'DELP_DRY', 'PEDGE_DRY', 'QI', 'QL', 'CMFMC', 'DTRAIN',
        'F_OF_PBL', 'DQRCU', 'DQRLSAN', 'PFICU', 'PFILSAN', 'PFLCU', 
        'PFLLSAN', 'TAUCLI', 'TAUCLW', 'AIRNUMDEN', 'AVGW', 'SPHU',
        'TV', 'DAIRMASS', 'F_UNDER_PBLTOP'
    ]
    
    # Define categorical 3D fields (different third dimension)
    categorical_3d_fields = [
        'SOILM',       # (nx,ny,nsoil) - soil layers
        'FRLANDUSE',   # (nx,ny,nlanduse) - land use categories
        'FRSOIL',      # (nx,ny,nsoil) - soil type fractions
        'FRLAI',       # (nx,ny,nlanduse) - LAI per land use category
        'FRZ0'         # (nx,ny,nlanduse) - roughness per land use category
    ]
    
    # Define the 2D fields we want to include
    wanted_2d_fields = [
        'T2M', 'PS', 'USTAR', 'U10M', 'V10M', 'PBLH', 'Z0', 'LAI', 'FRLAND',
        'FROCEAN', 'FRSEAICE', 'AREA_M2', 'SST', 'TSKIN', 'TS', 'QV2M', 'SLP',
        'PHIS', 'SUNCOS', 'SWGDN', 'HFLUX', 'EFLUX', 'PRECCON', 'PRECLSC',
        'PRECANV', 'TO3', 'TROPP', 'TropHt'
    ]
    
    with open(output_file, 'w') as f:
        f.write("! Generated VirtualMetType cleanup procedure\n")
        f.write("! This macro should be included in the virtual_met_cleanup subroutine\n\n")
        
        f.write("      ! Nullify 3D atmospheric field pointers (do not deallocate - they point to MetState data)\n")
        for name, rank, dims in fields:
            if rank == 3 and name in wanted_3d_fields:
                f.write(f"      this%{name} => null()\n")
        
        f.write("\n      ! Nullify 3D categorical field pointers (do not deallocate - they point to MetState data)\n")
        for name, rank, dims in fields:
            if rank == 3 and name in categorical_3d_fields:
                f.write(f"      this%{name} => null()\n")
        
        f.write("\n      ! Reset 2D scalar fields to default values\n")
        for name, rank, dims in fields:
            if rank == 2 and name in wanted_2d_fields:
                f.write(f"      this%{name} = 0.0_fp\n")
        
        f.write("\n")

def write_virtualcolumn_populate(fields, output_file):
    """
    Write a macro for populating VirtualColumn with all MetState fields (old approach).

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    
    # Define categorical 3D fields that are not vertical profiles
    categorical_3d_fields = {
        'SOILM',       # (nx,ny,nsoil)
        'FRLANDUSE',   # (nx,ny,nlanduse)  
        'FRSOIL',      # (nx,ny,nsoil)
        'FRLAI',       # (nx,ny,nlanduse)
        'FRZ0'         # (nx,ny,nlanduse)
    }
    
    with open(output_file, 'w') as f:
        f.write("! Generated macro for populating VirtualColumn with MetState fields\n")
        f.write("! This macro should be included in the populate_virtual_column subroutine\n\n")
        
        # Add 3D vertical profile fields only (exclude categorical 3D fields)
        for name, rank, dims in fields:
            if rank == 3 and name not in categorical_3d_fields:
                f.write(f"! Add {name} 3D vertical profile field\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   call virtual_col%add_met_field('{name}', this%met_state%{name}(grid_i, grid_j, :), rc)\n")
                f.write(f"   if (rc /= CC_SUCCESS) return\n")
                f.write(f"end if\n\n")
        
        # Add 2D surface fields
        for name, rank, dims in fields:
            if rank == 2:
                f.write(f"! Add {name} 2D surface field\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   allocate(temp_column(1))\n")
                f.write(f"   temp_column(1) = this%met_state%{name}(grid_i, grid_j)\n")
                f.write(f"   call virtual_col%add_met_field('{name}', temp_column, rc)\n")
                f.write(f"   deallocate(temp_column)\n")
                f.write(f"   if (rc /= CC_SUCCESS) return\n")
                f.write(f"end if\n\n")
        
        # Add scalar fields
        for name, rank, dims in fields:
            if rank == 0:
                f.write(f"! Add {name} scalar field\n")
                f.write(f"allocate(temp_column(1))\n")
                f.write(f"temp_column(1) = this%met_state%{name}\n")
                f.write(f"call virtual_col%add_met_field('{name}', temp_column, rc)\n")
                f.write(f"deallocate(temp_column)\n")
                f.write(f"if (rc /= CC_SUCCESS) return\n\n")
        
        f.write("! Note: Categorical 3D fields (SOILM, FRLANDUSE, FRSOIL, FRLAI, FRZ0)\n")
        f.write("! are excluded as they require specific category indices, not column extraction.\n")
        f.write("! Processes needing these should access them directly from MetState.\n")
    """
    Write a macro for populating VirtualColumn with all MetState fields.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    
    # Define categorical 3D fields that are not vertical profiles
    # These have 3rd dimensions like nsoil, nlanduse, etc. rather than nz
    categorical_3d_fields = {
        'SOILM',       # (nx,ny,nsoil)
        'FRLANDUSE',   # (nx,ny,nlanduse)  
        'FRSOIL',      # (nx,ny,nsoil)
        'FRLAI',       # (nx,ny,nlanduse)
        'FRZ0'         # (nx,ny,nlanduse)
    }
    
    with open(output_file, 'w') as f:
        f.write("! Generated macro for populating VirtualColumn with MetState fields\n")
        f.write("! This macro should be included in the populate_virtual_column subroutine\n\n")
        
        # Add 3D vertical profile fields only (exclude categorical 3D fields)
        for name, rank, dims in fields:
            if rank == 3 and name not in categorical_3d_fields:
                f.write(f"! Add {name} 3D vertical profile field\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   call virtual_col%add_met_field('{name}', this%met_state%{name}(grid_i, grid_j, :), rc)\n")
                f.write(f"   if (rc /= CC_SUCCESS) return\n")
                f.write(f"end if\n\n")
        
        # Add 2D surface fields
        for name, rank, dims in fields:
            if rank == 2:
                f.write(f"! Add {name} 2D surface field\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   allocate(temp_column(1))\n")
                f.write(f"   temp_column(1) = this%met_state%{name}(grid_i, grid_j)\n")
                f.write(f"   call virtual_col%add_met_field('{name}', temp_column, rc)\n")
                f.write(f"   deallocate(temp_column)\n")
                f.write(f"   if (rc /= CC_SUCCESS) return\n")
                f.write(f"end if\n\n")
        
        # Add scalar fields
        for name, rank, dims in fields:
            if rank == 0:
                f.write(f"! Add {name} scalar field\n")
                f.write(f"allocate(temp_column(1))\n")
                f.write(f"temp_column(1) = this%met_state%{name}\n")
                f.write(f"call virtual_col%add_met_field('{name}', temp_column, rc)\n")
                f.write(f"deallocate(temp_column)\n")
                f.write(f"if (rc /= CC_SUCCESS) return\n\n")
        
        # Add note about excluded categorical fields
        f.write("! Note: Categorical 3D fields (SOILM, FRLANDUSE, FRSOIL, FRLAI, FRZ0)\n")
        f.write("! are excluded as they require specific category indices, not column extraction.\n")
        f.write("! Processes needing these should access them directly from MetState.\n")

def write_scalar_accessor(fields, output_file):
    """
    Write a macro for scalar (rank-0) field accessors.

    Parameters
    ----------
    fields : list of (name, rank, dims)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    # Only scalar fields (rank 0)
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 0:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   scalar_val = this%{name}\n")
        f.write("case default\n")
        f.write("   scalar_val = 0.0_fp\n")
    """
    Write a macro for scalar (rank-0) field accessors.

    Parameters
    ----------
    fields : list of (name, rank, dims)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    # Only scalar fields (rank 0)
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 0:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   scalar_val = this%{name}\n")
        f.write("case default\n")
        f.write("   scalar_val = 0.0_fp\n")

def main():
    """
    Main entry point for the macro generator script.

    Parses command-line arguments and dispatches to the appropriate macro writer.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate Fortran macros for MetStateType fields.",
        epilog="Example: generate_metstate_macros.py allocate metstate_mod.F90 metstate_allocate_fields.inc"
    )
    parser.add_argument('mode', choices=['accessor', 'allocate', 'deallocate', 'column_accessor', '2d_scalar_accessor', 'scalar_accessor', 'virtualcolumn_populate', 'virtualmet_populate', 'virtualmet_type', 'virtualmet_cleanup'],
                        help="Type of macro to generate.")
    parser.add_argument('input_f90', help="Input Fortran file containing MetStateType.")
    parser.add_argument('output_inc', help="Output .inc file to write.")
    args = parser.parse_args()

    fields = parse_metstate_type(args.input_f90)
    if not fields:
        print(f"Warning: No fields found in {args.input_f90}.")

    if args.mode == 'accessor':
        write_accessor(fields, args.output_inc)
    elif args.mode == 'allocate':
        write_allocate(fields, args.output_inc)
    elif args.mode == 'deallocate':
        write_deallocate(fields, args.output_inc)
    elif args.mode == 'column_accessor':
        write_column_accessor(fields, args.output_inc)
    elif args.mode == '2d_scalar_accessor':
        write_2d_scalar_accessor(fields, args.output_inc)
    elif args.mode == 'scalar_accessor':
        write_scalar_accessor(fields, args.output_inc)
    elif args.mode == 'virtualcolumn_populate':
        write_virtualcolumn_populate(fields, args.output_inc)
    elif args.mode == 'virtualmet_populate':
        write_virtualmet_populate(fields, args.output_inc)
    elif args.mode == 'virtualmet_type':
        write_virtualmet_type(fields, args.output_inc)
    elif args.mode == 'virtualmet_cleanup':
        write_virtualmet_cleanup(fields, args.output_inc)
    else:
        parser.error(f"Unknown mode: {args.mode}")

if __name__ == "__main__":
    main()
