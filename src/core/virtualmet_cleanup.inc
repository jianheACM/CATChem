! Generated VirtualMetType cleanup procedure
! This macro should be included in the virtual_met_cleanup subroutine

      ! Nullify 3D atmospheric field pointers (do not deallocate - they point to MetState data)
      this%F_OF_PBL => null()
      this%F_UNDER_PBLTOP => null()
      this%CLDF => null()
      this%CMFMC => null()
      this%DQRCU => null()
      this%DQRLSAN => null()
      this%DTRAIN => null()
      this%QI => null()
      this%QL => null()
      this%PFICU => null()
      this%PFILSAN => null()
      this%PFLCU => null()
      this%PFLLSAN => null()
      this%TAUCLI => null()
      this%TAUCLW => null()
      this%Z => null()
      this%ZMID => null()
      this%BXHEIGHT => null()
      this%QV => null()
      this%T => null()
      this%THETA => null()
      this%TV => null()
      this%V => null()
      this%U => null()
      this%OMEGA => null()
      this%RH => null()
      this%SPHU => null()
      this%AIRDEN => null()
      this%AIRNUMDEN => null()
      this%MAIRDEN => null()
      this%AVGW => null()
      this%DELP => null()
      this%DELP_DRY => null()
      this%DAIRMASS => null()
      this%AIRVOL => null()
      this%PEDGE_DRY => null()
      this%PEDGE => null()
      this%PMID => null()
      this%PMID_DRY => null()

      ! Nullify 3D categorical field pointers (do not deallocate - they point to MetState data)
      this%SOILM => null()
      this%FRLANDUSE => null()
      this%FRSOIL => null()
      this%FRLAI => null()
      this%FRZ0 => null()

      ! Reset 2D scalar fields to default values
      this%AREA_M2 = 0.0_fp
      this%FRLAND = 0.0_fp
      this%FROCEAN = 0.0_fp
      this%FRSEAICE = 0.0_fp
      this%LAI = 0.0_fp
      this%SUNCOS = 0.0_fp
      this%SWGDN = 0.0_fp
      this%EFLUX = 0.0_fp
      this%HFLUX = 0.0_fp
      this%U10M = 0.0_fp
      this%USTAR = 0.0_fp
      this%V10M = 0.0_fp
      this%Z0 = 0.0_fp
      this%PBLH = 0.0_fp
      this%PRECANV = 0.0_fp
      this%PRECCON = 0.0_fp
      this%PRECLSC = 0.0_fp
      this%PHIS = 0.0_fp
      this%QV2M = 0.0_fp
      this%T2M = 0.0_fp
      this%TS = 0.0_fp
      this%TSKIN = 0.0_fp
      this%SST = 0.0_fp
      this%SLP = 0.0_fp
      this%PS = 0.0_fp
      this%TO3 = 0.0_fp
      this%TROPP = 0.0_fp
      this%TropHt = 0.0_fp

