# CATChem Process Generator Demo

This is a hands-on demonstration of the CATChem process generator.

## Demo: Creating a Particle Settling Process

Let's create a complete particle settling process with two schemes.

### Step 1: Run the Generator

```bash
cd /path/to/catchem
python util/catchem_generate_process.py \
    --name Settling \
    --type transformation \
    --schemes stokes,intermediate_reynolds \
    --species PM25,PM10,DUST_1,DUST_2,DUST_3 \
    --description "Gravitational settling of atmospheric particles"
```

### Step 2: What Gets Generated

The generator creates the following files:

```
src/process/settling/
├── SettlingProcess_Mod.F90                    # Main process module
├── CMakeLists.txt                             # Build configuration
└── schemes/
    ├── SettlingStokesScheme_Mod.F90           # Stokes scheme
    └── SettlingIntermediateReynoldsScheme_Mod.F90  # Intermediate Reynolds scheme

tests/
└── test_settling.F90                          # Unit tests

docs/processes/settling/
├── index.md                                   # Main documentation
├── schemes.md                                 # Scheme descriptions
└── examples.md                                # Usage examples
```

### Step 3: Generated Code Structure

**Main Process Module (`SettlingProcess_Mod.F90`):**
```fortran
module SettlingProcess_Mod
  use ProcessInterface_Mod
  use StateContainer_Mod
  use ErrorManager_Mod
  implicit none
  private

  type, extends(ProcessInterface), public :: SettlingProcessType
    private
    character(len=32) :: selected_scheme = 'stokes'
    ! Process-specific members will be added here
  contains
    procedure, public :: init => settling_init
    procedure, public :: run => settling_run
    procedure, public :: finalize => settling_finalize
    procedure, private :: validate_chem_inputs
    procedure, private :: update_chemistry
    procedure, private :: check_mass_conservation
  end type SettlingProcessType

contains

  subroutine settling_init(this, container, rc)
    class(SettlingProcessType), intent(inout) :: this
    type(StateContainerType), intent(inout) :: container
    integer, intent(out) :: rc

    ! Implementation will be added here
    rc = 0
  end subroutine settling_init

  subroutine settling_run(this, container, rc)
    class(SettlingProcessType), intent(inout) :: this
    type(StateContainerType), intent(inout) :: container
    integer, intent(out) :: rc

    ! Implementation will be added here
    rc = 0
  end subroutine settling_run

  subroutine settling_finalize(this, rc)
    class(SettlingProcessType), intent(inout) :: this
    integer, intent(out) :: rc

    ! Implementation will be added here
    rc = 0
  end subroutine settling_finalize

end module SettlingProcess_Mod
```

**Stokes Scheme (`SettlingStokesScheme_Mod.F90`):**
```fortran
module SettlingStokesScheme_Mod
  use precision_mod
  implicit none
  private

  public :: stokes_calculate

contains

  subroutine stokes_calculate(container, process, rc)
    use StateContainer_Mod
    use SettlingProcess_Mod

    type(StateContainerType), intent(inout) :: container
    class(SettlingProcessType), intent(inout) :: process
    integer, intent(out) :: rc

    ! Stokes settling calculation implementation goes here
    ! v_settling = (2 * r^2 * rho_p * g) / (9 * mu)

    rc = 0
  end subroutine stokes_calculate

end module SettlingStokiesScheme_Mod
```

### Step 4: Configuration

The process can be configured via YAML:

```yaml
processes:
  - name: settling
    enabled: true
    scheme: stokes  # or intermediate_reynolds
    species: [PM25, PM10, DUST_1, DUST_2, DUST_3]
    parameters:
      particle_density: 2650.0  # kg/m³
      reference_diameter: 1.0e-6  # m
    diagnostics:
      - settling_velocity
      - particle_flux
```

### Step 5: Testing

Generated unit tests (`test_settling.F90`):
```fortran
program test_settling
  use SettlingProcess_Mod
  use StateContainer_Mod
  use testing_mod
  implicit none

  type(SettlingProcessType) :: process
  type(StateContainerType) :: container
  integer :: rc

  call test_settling_init()
  call test_settling_run()
  call test_settling_finalize()

contains

  subroutine test_settling_init()
    ! Test initialization
    call process%init(container, rc)
    call assert_equals(rc, 0, "Settling init should succeed")
  end subroutine

  ! Additional tests...

end program test_settling
```

### Step 6: Documentation

The generator creates comprehensive documentation at `docs/processes/settling/`:

**`index.md`** - Main process documentation with scientific background
**`schemes.md`** - Detailed scheme descriptions and parameters
**`examples.md`** - Configuration and usage examples

### Step 7: Next Steps

After generation, you would:

1. **Implement the physics** in the scheme modules
2. **Add real parameters** and validation logic
3. **Enhance tests** with specific test cases
4. **Update documentation** with scientific references

## Try It Yourself!

Run the generator with different parameters to see how it adapts:

```bash
# Simple emission process
python util/catchem_generate_process.py \
    --name BiogenicEmissions \
    --type emission \
    --schemes megan,beis \
    --species ISOP,TERP,NO

# Complex multi-phase chemistry
python util/catchem_generate_process.py \
    --name AqueousChemistry \
    --type multiphase_chemistry \
    --phases gas,liquid \
    --solver micm \
    --schemes cb6_aqueous

# Transport process
python util/catchem_generate_process.py \
    --name Advection \
    --type transport \
    --schemes upwind,ppm
```

The generator adapts the templates and structure based on the process type, creating appropriate patterns for each use case.
