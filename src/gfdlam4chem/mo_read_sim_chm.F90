      module MO_READ_SIM_CHM_MOD

use mo_errmsg, only : errmsg

implicit none
private
public :: READ_SIM_CHM

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      CONTAINS
        
      subroutine READ_SIM_CHM( sim_data_flsp, &
                               sim_file_cnt )
!--------------------------------------------------------
!            ... Initialize chemistry modules
!--------------------------------------------------------

      use CHEM_MODS_MOD,     only : explicit, implicit, rodas, grpcnt, &
                                nadv_mass, adv_mass, pcnstm1, &
                                drydep_cnt, drydep_lst, &
                                srfems_cnt, srfems_lst, &
                                hetcnt, het_lst, extcnt, extfrc_lst, &
                                rxt_alias_cnt, rxt_alias_lst, rxt_alias_map, &
                                ngrp, grp_mem_cnt, grp_lst
      use M_TRACNAME_MOD,    only : tracnam, natsnam

      implicit none

!--------------------------------------------------------
!            ... Dummy args
!--------------------------------------------------------
      integer, intent(out) :: sim_file_cnt
      character(len=32), intent(in) :: sim_data_flsp

!--------------------------------------------------------
!            ... Local variables
!--------------------------------------------------------
      integer, parameter :: inst = 1, avrg = 2
      integer, parameter :: max_hst_ind = 17
      integer  ::  ios, funit
      integer  ::  moz_file_cnt, match_file_cnt
      character(len=128) :: msg

!     funit = NAVU()
!--------------------------------------------------------
!            ... Open chem input unit
!--------------------------------------------------------
     funit = 66
     OPEN( unit = funit, &
           file = TRIM( sim_data_flsp ), &
           action='read', &
           status = 'old', &
           recl   = 2048, &
           iostat = ios )
     if( ios /= 0 ) then
        write(msg,*) ' READ_SIM_CHM: Failed to open file ',TRIM( sim_data_flsp )
        call ENDRUN(msg)
     end if
      
      if( explicit%clscnt > 0 ) then
         read(funit,'(4i4)',iostat=ios) explicit%cls_rxt_cnt(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read explicit cls_rxt_cnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) explicit%clsmap(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read explicit clscnt; error = ', ios
            call ENDRUN(msg)
         end if
      end if

      if( implicit%clscnt > 0 ) then
         read(funit,'(4i4)',iostat=ios) implicit%cls_rxt_cnt(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit cls_rxt_cnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) implicit%clsmap(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit clscnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) implicit%permute(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit permute; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) implicit%diag_map(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit diag_map; error = ', ios
            call ENDRUN(msg)
         end if
      end if

     if( rodas%clscnt > 0 ) then
         read(funit,'(4i4)',iostat=ios) rodas%cls_rxt_cnt(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas cls_rxt_cnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rodas%clsmap(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas clscnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rodas%permute(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas permute; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rodas%diag_map(:)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas diag_map; error = ', ios
            call ENDRUN(msg)
         end if
      end if

      if( pcnstm1 > 0 ) then
         read(funit,*,iostat=ios) adv_mass(:pcnstm1)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read adv_mass; error = ', ios
            call ENDRUN(msg)
         end if
      end if

      if( grpcnt > 0 ) then
         read(funit,*,iostat=ios) nadv_mass(:grpcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read nadv_mass; error = ', ios
            call ENDRUN(msg)
         end if
      end if

      if( pcnstm1 > 0 ) then
         read(funit,'(10a8)',iostat=ios) tracnam(:pcnstm1)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read tracnam; error = ', ios
            call ENDRUN(msg)
         end if
      end if

      if( grpcnt > 0 ) then
         read(funit,'(i4)',iostat=ios) ngrp
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read ngrp; error = ',ios
            call ENDRUN(msg)
         end if
         allocate( grp_mem_cnt(ngrp),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate grp_mem_cnt; error = ',ios
            call ENDRUN(msg)
         end if
         allocate( grp_lst(ngrp),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate grp_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) grp_mem_cnt(:ngrp)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read grp_mem_cnt; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) grp_lst(:ngrp)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read grp_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) natsnam(1:grpcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read natsnam; error = ',ios
            call ENDRUN(msg)
         end if
      end if

      read(funit,'(i4)',iostat=ios) srfems_cnt
      if( ios /= 0 ) then
         write(msg,*) 'READ_SIM_CHM: Failed to read srfems_cnt; error = ',ios
            call ENDRUN(msg)
      end if
      if( srfems_cnt > 0 ) then
         allocate( srfems_lst(srfems_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate srfems_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) srfems_lst(1:srfems_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read srfems_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if

      read(funit,'(i4)',iostat=ios) drydep_cnt
      if( ios /= 0 ) then
         write(msg,*) 'READ_SIM_CHM: Failed to read drydep_cnt; error = ',ios
            call ENDRUN(msg)
      end if
      if( drydep_cnt > 0 ) then
         allocate( drydep_lst(drydep_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate drydep_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) drydep_lst(1:drydep_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read drydep_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if

      if( hetcnt > 0 ) then
         read(funit,'(10a8)',iostat=ios) het_lst(1:hetcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read het_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if

      if( extcnt > 0 ) then
         read(funit,'(10a8)',iostat=ios) extfrc_lst(1:extcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read extfrc_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if

      read(funit,'(i4)',iostat=ios) rxt_alias_cnt
      if( ios /= 0 ) then
         write(msg,*) 'READ_SIM_CHM: Failed to read rxt_alias_cnt; error = ',ios
            call ENDRUN(msg)
      end if
      if( rxt_alias_cnt > 0 ) then
         allocate( rxt_alias_lst(rxt_alias_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate rxt_alias_lst; error = ',ios
            call ENDRUN(msg)
         end if
         allocate( rxt_alias_map(rxt_alias_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate rxt_alias_map; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(5a16)',iostat=ios) rxt_alias_lst(1:rxt_alias_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rxt_alias_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rxt_alias_map(1:rxt_alias_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rxt_alias_map; error = ',ios
            call ENDRUN(msg)
         end if
      end if

     close(funit)
     ! call close_file(funit,dist=.true.)

!      write(*,*) '---------------------------------------------------------------------------------'
!      write(*,*) ' '

      end subroutine READ_SIM_CHM

      subroutine ENDRUN(msg)
         character(len=128), intent(in) :: msg
         call errmsg("Error in mo_read_sim_chm",msg,.true.)
      end subroutine ENDRUN        

      end module MO_READ_SIM_CHM_MOD
