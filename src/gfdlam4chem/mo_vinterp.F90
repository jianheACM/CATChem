module  mo_vinterp

  use catchem_constants, only : kind_chem
  use mo_errmsg,             only : errmsg

  implicit none

  public

contains

!---------------------------------------------------------------------
!> \brief interp_weighted_scalar_2D receives the variables grdin,
!!            grdout, and datin as inputs and returns datout.
!!
!! \param [in] <grdin> No description
!! \param [in] <grdout> No description
!! \param [in] <datin> No description
!! \param [out] <datout> No description
 subroutine interp_weighted_scalar_2D (grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout
real, intent(in),  dimension(:,:) :: datin
real, intent(out), dimension(:,:) :: datout

integer :: j, k, n

if (size(grdin(:)).ne. (size(datin,1)+1)) &
 call errmsg('interp_weighted_scalar: ', 'input data and pressure do not have the same number of levels', .true.)
if (size(grdout(:)).ne. (size(datout,1 )+1)) &
 call errmsg('interp_weighted_scalar: ', 'output data and pressure do not have the same number of levels',.true.)

  do k = 1, size(datout,1 )
   datout(k,:) = 0.0

     do j = 1, size(datin,1 )

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdin(j+1)-grdout(k))
          end do

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdout(k+1)-grdin(j))
          end do

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdin(j+1)-grdin(j))
          end do

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

          do n= 1, size(datin,2)
          datout(k,n) = datout(k,n) + datin(j,n)*(grdout(k+1)-grdout(k))

          end do
        endif

     enddo

     do n= 1, size(datin,2)
       datout(k,n) = datout(k,n)/(grdout(k+1)-grdout(k))
     end do

  enddo

end subroutine interp_weighted_scalar_2D

!
!---------------------------------------------------------------------
!> \brief interp_weighted_scalar_1D receives the variables grdin,
!!        grdout, and datin as inputs and returns datout.
!!
!! \param [in] <grdin> No description
!! \param [in] <grdout> No description
!! \param [in] <datin> No description
!! \param [out] <datout> No description
 subroutine interp_weighted_scalar_1D (grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout, datin
real, intent(out), dimension(:) :: datout

integer :: j, k

if (size(grdin(:)).ne. (size(datin(:))+1)) &
 call errmsg('interp_weighted_scalar: ', 'input data and pressure do not have the same number of levels',.true.)
if (size(grdout(:)).ne. (size(datout(:))+1)) &
 call errmsg('interp_weighted_scalar: ', 'output data and pressure do not have the same number of levels',.true.)

  do k = 1, size(datout(:))
   datout(k) = 0.0

     do j = 1, size(datin(:))

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdout(k))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdin(j))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdin(j))

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdout(k))

        endif

     enddo

     datout(k) = datout(k)/(grdout(k+1)-grdout(k))

  enddo

end subroutine interp_weighted_scalar_1D

!
!#################################################################
!
!---------------------------------------------------------------------
!> \brief interp_linear receives the variables grdin,
!!            grdout, and datin as inputs and returns a linear
!!            interpolation.
!!
!! \param [in] <grdin> No description
!! \param [in] <grdout> No description
!! \param [in] <datin> No description
!! \param [out] <datout> No description
subroutine interp_linear ( grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout, datin
real, intent(out), dimension(:) :: datout

integer :: j, k, n
real    :: wt


if (size(grdin(:)).ne. (size(datin(:))+1)) &
 call errmsg('interp_linear: ', 'input data and pressure do not have the same number of levels',.true.)
if (size(grdout(:)).ne. (size(datout(:))+1)) &
 call errmsg('interp_linear: ', 'output data and pressure do not have the same number of levels',.true.)


  n = size(grdin(:))

  do k= 1, size(datout(:))

   ! ascending grid values
     if (grdin(1) < grdin(n)) then
         do j = 2, size(grdin(:))-1
           if (grdout(k) <= grdin(j)) exit
         enddo
   ! descending grid values
     else
         do j = size(grdin(:)), 3, -1
           if (grdout(k) <= grdin(j-1)) exit
         enddo
     endif

   ! linear interpolation
     wt = (grdout(k)-grdin(j-1)) / (grdin(j)-grdin(j-1))
!print '(a,2i3,4f6.1)', 'k,j=',k,j,grdout(k),grdin(j-1),grdin(j),wt
   ! constant value extrapolation
   ! wt = min(max(wt,0.),1.)

     datout(k) = (1.-wt)*datin(j-1) + wt*datin(j)

  enddo

end subroutine interp_linear
!


  subroutine vinterp(chem_in,nspecies,nlv, p_in, pphy_mb, chem_out)

    real(kind=kind_chem), dimension(1:nlv,1:nspecies) :: chem_in
    real(kind=kind_chem), dimension(1:nlv), intent(in) :: p_in
    real(kind=kind_chem), intent(in) :: pphy_mb
    integer, intent(in) :: nspecies, nlv 
    real (kind=kind_chem), dimension(1:nspecies),intent(out) :: chem_out

    integer :: l, ll, n
    real (kind=kind_chem) :: pu, pl, pwant, aln

    chem_out = 0.

    do ll = 2, nlv
      l = ll
      if (p_in(l) < pphy_mb) exit
    enddo

    pu=alog(p_in(l))
    pl=alog(p_in(l-1))
    pwant=alog(pphy_mb)

    do n = 1, nspecies
      if (pwant > pl) then
        chem_out(n) = chem_in(l,n)
      else
        aln=(chem_in(l,n)*(pwant-pl)+            &
           chem_in(l-1,n)*(pu-pwant))/(pu-pl)
        chem_out(n) = aln
      endif
    enddo
  end subroutine vinterp

end module mo_vinterp
