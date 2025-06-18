!> \file emisstate_mod.F90
!! \brief Module for emission state variables
!!
!! This module contains the EmisStateType derived type and related procedures
!! for managing emission data in CATChem. It follows the modern state container
!! pattern with proper initialization, cleanup, and validation procedures.
!!
!! \ingroup core_modules
!!!>
MODULE EmisState_Mod

   !=========================================================================
   ! Module uses
   !=========================================================================
   USE Error_Mod
   USE Precision_Mod
   USE species_mod, only: SpeciesType

   IMPLICIT NONE
   PRIVATE

   !=========================================================================
   ! Public interfaces
   !=========================================================================
   PUBLIC :: EmisStateType
   PUBLIC :: EmisSpeciesType
   PUBLIC :: EmisCategoryType
   ! All legacy functions removed - use modern type-bound procedures instead:
   ! - EmisState%init() instead of legacy allocation
   ! - EmisState%cleanup() instead of EmisState_CleanUp()
   ! - EmisState%find_species() for species mapping
   ! - EmisState%apply_emissions() for applying emissions

   !=========================================================================
   ! Derived types for emission species
   !=========================================================================

   !> \brief Derived type for individual emission species
   !!
   !! Contains all data for a single emitted species including fluxes,
   !! plume rise parameters, and mapping information.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: EmisSpeciesType
      CHARACTER(LEN=31)             :: name             ! Species name
      INTEGER                       :: nEmisMap         ! Number of emission maps
      INTEGER,          ALLOCATABLE :: EmisMapIndex(:)  ! Emission map indices
      REAL(fp),         ALLOCATABLE :: Flux(:)          ! Emission flux [kg/m2/s]
      INTEGER                       :: plumerise        ! Plume rise option
      INTEGER                       :: nPlmSrc          ! Number of plume sources
      REAL(fp),         ALLOCATABLE :: frp(:)           ! Fire radiative power [MW]
      REAL(fp),         ALLOCATABLE :: PlmSrcFlx(:)     ! Plume source flux [kg/s]
      REAL(fp),         ALLOCATABLE :: PlmRiseHgt(:)    ! Plume rise height [m]
      REAL(fp),         ALLOCATABLE :: STKDM(:)         ! Stack diameter [m]
      REAL(fp),         ALLOCATABLE :: STKHT(:)         ! Stack height [m]
      REAL(fp),         ALLOCATABLE :: STKTK(:)         ! Stack temperature [K]
      REAL(fp),         ALLOCATABLE :: STKVE(:)         ! Stack velocity [m/s]

   CONTAINS
      PROCEDURE :: init => emisspecies_init
      PROCEDURE :: cleanup => emisspecies_cleanup
      PROCEDURE :: validate => emisspecies_validate

   END TYPE EmisSpeciesType

   !> \brief Derived type for emission categories
   !!
   !! Contains data for a single emission category including all species
   !! within that category and category-specific parameters.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: EmisCategoryType
      CHARACTER(LEN=31)                   :: Name           ! Category name
      INTEGER                             :: nSpecies       ! Number of species
      INTEGER                             :: nPlumerise     ! Number of plume rise species
      TYPE(EmisSpeciesType), ALLOCATABLE  :: Species(:)     ! Species array

   CONTAINS
      PROCEDURE :: init => emiscategory_init
      PROCEDURE :: cleanup => emiscategory_cleanup
      PROCEDURE :: validate => emiscategory_validate

   END TYPE EmisCategoryType

   !> \brief Derived type for emission state
   !!
   !! Main container for all emission data including categories, species,
   !! and global emission parameters.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: EmisStateType
      CHARACTER(LEN=4)                    :: State              = 'EMIS'
      INTEGER                             :: nCats               ! Number of categories
      INTEGER                             :: nEmisTotalPlumerise ! Total plume rise emissions
      TYPE(EmisCategoryType), ALLOCATABLE :: Cats(:)            ! Category array

   CONTAINS
      PROCEDURE :: init => emisstate_init
      PROCEDURE :: cleanup => emisstate_cleanup_method
      PROCEDURE :: validate => emisstate_validate

   END TYPE EmisStateType

CONTAINS

   !=========================================================================
   ! EmisSpeciesType procedures
   !=========================================================================

   !> \brief Initialize an emission species
   !!
   !! \param[inout] this The EmisSpeciesType object to initialize
   !! \param[in] species_name Name of the species
   !! \param[in] num_levels Number of vertical levels
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   !! \param[in] num_plume_sources Number of plume sources (optional)
   !!!>
   SUBROUTINE emisspecies_init(this, species_name, num_levels, error_mgr, rc, num_plume_sources)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_MEMORY_ALLOCATION, ERROR_INVALID_INPUT

      CLASS(EmisSpeciesType), INTENT(INOUT) :: this
      CHARACTER(LEN=*),       INTENT(IN)    :: species_name
      INTEGER,                INTENT(IN)    :: num_levels
      TYPE(ErrorManagerType), POINTER, INTENT(INOUT) :: error_mgr
      INTEGER,                INTENT(OUT)   :: rc
      INTEGER, OPTIONAL,      INTENT(IN)    :: num_plume_sources

      INTEGER :: nPlmSrc

      ! Initialize
      call error_mgr%push_context('emisspecies_init', 'emisstate_mod.F90')
      rc = CC_SUCCESS

      ! Validate input
      if (num_levels <= 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Number of levels must be positive', rc)
         call error_mgr%pop_context()
         return
      endif

      ! Set defaults
      nPlmSrc = 0
      if (present(num_plume_sources)) nPlmSrc = num_plume_sources

      ! Set basic properties
      this%name = TRIM(species_name)
      this%nEmisMap = 0
      this%plumerise = 0
      this%nPlmSrc = nPlmSrc

      ! Allocate flux array
      if (allocated(this%Flux)) deallocate(this%Flux)
      allocate(this%Flux(num_levels), stat=rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                     'Error allocating Flux array for species ' // TRIM(species_name), rc)
         call error_mgr%pop_context()
         return
      endif
      this%Flux = 0.0_fp

      ! Allocate plume rise arrays if needed
      if (nPlmSrc > 0) then
         if (allocated(this%frp)) deallocate(this%frp)
         if (allocated(this%PlmSrcFlx)) deallocate(this%PlmSrcFlx)
         if (allocated(this%PlmRiseHgt)) deallocate(this%PlmRiseHgt)
         if (allocated(this%STKDM)) deallocate(this%STKDM)
         if (allocated(this%STKHT)) deallocate(this%STKHT)
         if (allocated(this%STKTK)) deallocate(this%STKTK)
         if (allocated(this%STKVE)) deallocate(this%STKVE)

         allocate(this%frp(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating frp array', rc)
            call error_mgr%pop_context()
            return
         endif

         allocate(this%PlmSrcFlx(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating PlmSrcFlx array', rc)
            call error_mgr%pop_context()
            return
         endif

         allocate(this%PlmRiseHgt(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating PlmRiseHgt array', rc)
            call error_mgr%pop_context()
            return
         endif

         allocate(this%STKDM(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating STKDM array', rc)
            call error_mgr%pop_context()
            return
         endif

         allocate(this%STKHT(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating STKHT array', rc)
            call error_mgr%pop_context()
            return
         endif

         allocate(this%STKTK(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating STKTK array', rc)
            call error_mgr%pop_context()
            return
         endif

         allocate(this%STKVE(nPlmSrc), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating STKVE array', rc)
            call error_mgr%pop_context()
            return
         endif

         ! Initialize arrays
         this%frp = 0.0_fp
         this%PlmSrcFlx = 0.0_fp
         this%PlmRiseHgt = 0.0_fp
         this%STKDM = 0.0_fp
         this%STKHT = 0.0_fp
         this%STKTK = 0.0_fp
         this%STKVE = 0.0_fp
      endif

      call error_mgr%pop_context()

   END SUBROUTINE emisspecies_init

   !> \brief Clean up an emission species
   !!
   !! \param[inout] this The EmisSpeciesType object to clean up
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emisspecies_cleanup(this, rc)
      CLASS(EmisSpeciesType), INTENT(INOUT) :: this
      INTEGER,                INTENT(OUT)   :: rc

      rc = CC_SUCCESS

      if (allocated(this%EmisMapIndex)) deallocate(this%EmisMapIndex)
      if (allocated(this%Flux)) deallocate(this%Flux)
      if (allocated(this%frp)) deallocate(this%frp)
      if (allocated(this%PlmSrcFlx)) deallocate(this%PlmSrcFlx)
      if (allocated(this%PlmRiseHgt)) deallocate(this%PlmRiseHgt)
      if (allocated(this%STKDM)) deallocate(this%STKDM)
      if (allocated(this%STKHT)) deallocate(this%STKHT)
      if (allocated(this%STKTK)) deallocate(this%STKTK)
      if (allocated(this%STKVE)) deallocate(this%STKVE)

      this%name = ''
      this%nEmisMap = 0
      this%plumerise = 0
      this%nPlmSrc = 0

   END SUBROUTINE emisspecies_cleanup

   !> \brief Validate an emission species
   !!
   !! \param[in] this The EmisSpeciesType object to validate
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emisspecies_validate(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_INVALID_STATE

      CLASS(EmisSpeciesType), INTENT(IN)  :: this
      TYPE(ErrorManagerType), POINTER, INTENT(INOUT) :: error_mgr
      INTEGER,                INTENT(OUT) :: rc

      call error_mgr%push_context('emisspecies_validate', 'emisstate_mod.F90')
      rc = CC_SUCCESS

      if (LEN_TRIM(this%name) == 0) then
         call error_mgr%report_error(ERROR_INVALID_STATE, &
                                     'Species name is empty', rc)
         call error_mgr%pop_context()
         return
      endif

      if (.not. allocated(this%Flux)) then
         call error_mgr%report_error(ERROR_INVALID_STATE, &
                                     'Flux array not allocated for species ' // TRIM(this%name), rc)
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()

   END SUBROUTINE emisspecies_validate

   !=========================================================================
   ! EmisCategoryType procedures
   !=========================================================================

   !> \brief Initialize an emission category
   !!
   !! \param[inout] this The EmisCategoryType object to initialize
   !! \param[in] category_name Name of the category
   !! \param[in] num_species Number of species in this category
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emiscategory_init(this, category_name, num_species, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_MEMORY_ALLOCATION, ERROR_INVALID_INPUT

      CLASS(EmisCategoryType), INTENT(INOUT) :: this
      CHARACTER(LEN=*),        INTENT(IN)    :: category_name
      INTEGER,                 INTENT(IN)    :: num_species
      TYPE(ErrorManagerType), POINTER, INTENT(INOUT) :: error_mgr
      INTEGER,                 INTENT(OUT)   :: rc

      ! Initialize
      call error_mgr%push_context('emiscategory_init', 'emisstate_mod.F90')
      rc = CC_SUCCESS

      if (num_species <= 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Number of species must be positive', rc)
         call error_mgr%pop_context()
         return
      endif

      ! Set basic properties
      this%Name = TRIM(category_name)
      this%nSpecies = num_species
      this%nPlumerise = 0

      ! Allocate species array
      if (allocated(this%Species)) deallocate(this%Species)
      if (num_species > 0) then
         allocate(this%Species(num_species), stat=rc)
         if (rc /= CC_SUCCESS) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                        'Error allocating Species array for category ' // TRIM(category_name), rc)
            call error_mgr%pop_context()
            return
         endif
      endif

      call error_mgr%pop_context()

   END SUBROUTINE emiscategory_init

   !> \brief Clean up an emission category
   !!
   !! \param[inout] this The EmisCategoryType object to clean up
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emiscategory_cleanup(this, rc)
      CLASS(EmisCategoryType), INTENT(INOUT) :: this
      INTEGER,                 INTENT(OUT)   :: rc

      INTEGER :: i

      rc = CC_SUCCESS

      ! Clean up all species
      if (allocated(this%Species)) then
         do i = 1, size(this%Species)
            call this%Species(i)%cleanup(rc)
            if (rc /= CC_SUCCESS) return
         end do
         deallocate(this%Species)
      endif

      this%Name = ''
      this%nSpecies = 0
      this%nPlumerise = 0

   END SUBROUTINE emiscategory_cleanup

   !> \brief Validate an emission category
   !!
   !! \param[in] this The EmisCategoryType object to validate
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emiscategory_validate(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_INVALID_STATE

      CLASS(EmisCategoryType), INTENT(IN)  :: this
      TYPE(ErrorManagerType), POINTER, INTENT(INOUT) :: error_mgr
      INTEGER,                 INTENT(OUT) :: rc

      INTEGER :: i

      call error_mgr%push_context('emiscategory_validate', 'emisstate_mod.F90')
      rc = CC_SUCCESS

      if (LEN_TRIM(this%Name) == 0) then
         call error_mgr%report_error(ERROR_INVALID_STATE, &
                                     'Category name is empty', rc)
         call error_mgr%pop_context()
         return
      endif

      if (this%nSpecies > 0) then
         if (.not. allocated(this%Species)) then
            call error_mgr%report_error(ERROR_INVALID_STATE, &
                                        'Species array not allocated for category ' // TRIM(this%Name), rc)
            call error_mgr%pop_context()
            return
         endif

         ! Validate all species
         do i = 1, this%nSpecies
            call this%Species(i)%validate(error_mgr, rc)
            if (rc /= CC_SUCCESS) return
         end do
      endif

   END SUBROUTINE emiscategory_validate

   !=========================================================================
   ! EmisStateType procedures
   !=========================================================================

   !> \brief Initialize an emission state
   !!
   !! \param[inout] this The EmisStateType object to initialize
   !! \param[in] num_categories Number of emission categories
   !! \param[in] error_mgr Error manager for handling errors
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emisstate_init(this, num_categories, error_mgr, rc)
      CLASS(EmisStateType), INTENT(INOUT) :: this
      INTEGER,              INTENT(IN)    :: num_categories
      TYPE(ErrorManagerType), POINTER     :: error_mgr
      INTEGER,              INTENT(OUT)   :: rc

      INTEGER :: alloc_status

      ! Initialize
      rc = CC_SUCCESS

      ! Set basic properties
      this%nCats = num_categories
      this%nEmisTotalPlumerise = 0

      ! Allocate categories array
      if (allocated(this%Cats)) deallocate(this%Cats)
      if (num_categories > 0) then
         allocate(this%Cats(num_categories), stat=alloc_status)
         if (alloc_status /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                      'Error allocating Cats array', &
                                      rc, 'emisstate_init')
            rc = CC_FAILURE
            return
         endif
      endif

   END SUBROUTINE emisstate_init

   !> \brief Clean up an emission state
   !!
   !! \param[inout] this The EmisStateType object to clean up
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emisstate_cleanup_method(this, rc)
      CLASS(EmisStateType), INTENT(INOUT) :: this
      INTEGER,              INTENT(OUT)   :: rc

      INTEGER :: i

      rc = CC_SUCCESS

      ! Clean up all categories
      if (allocated(this%Cats)) then
         do i = 1, size(this%Cats)
            call this%Cats(i)%cleanup(rc)
            if (rc /= CC_SUCCESS) return
         end do
         deallocate(this%Cats)
      endif

      this%nCats = 0
      this%nEmisTotalPlumerise = 0

   END SUBROUTINE emisstate_cleanup_method

   !> \brief Validate an emission state
   !!
   !! \param[in] this The EmisStateType object to validate
   !! \param[in] error_mgr Error manager for handling errors
   !! \param[out] rc Return code
   !!!>
   SUBROUTINE emisstate_validate(this, error_mgr, rc)
      CLASS(EmisStateType), INTENT(IN)  :: this
      TYPE(ErrorManagerType), POINTER   :: error_mgr
      INTEGER,              INTENT(OUT) :: rc

      INTEGER :: i

      rc = CC_SUCCESS

      if (this%nCats > 0) then
         if (.not. allocated(this%Cats)) then
            call error_mgr%report_error(ERROR_INVALID_STATE, &
                                      'Cats array not allocated', &
                                      rc, 'emisstate_validate')
            rc = CC_FAILURE
            return
         endif

         ! Validate all categories
         do i = 1, this%nCats
            call this%Cats(i)%validate(error_mgr, rc)
            if (rc /= CC_SUCCESS) return
         end do
      endif

   END SUBROUTINE emisstate_validate

END MODULE EmisState_Mod
