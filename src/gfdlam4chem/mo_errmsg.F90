module  mo_errmsg

  implicit none

  public

contains

  subroutine errmsg(routine, messag, fatal)
!
! print out a warning or error message;  abort if error
!

      implicit none
      logical       fatal
      logical, save :: once
      data once / .false. /
      character*1000 msg
      character*(*) routine
      character*(*) messag
      integer, save :: maxmsg, nummsg
      data nummsg / 0 /,  maxmsg / 1000 /
!
!
      if ( fatal )  then
             write( msg, '(a)' )   &
                  routine // messag
      end if
!
      nummsg = nummsg + 1
      if ( nummsg.gt.maxmsg )  then
         if ( .not.once )then
            write( msg, '(a)' )   &
             routine // ' too many warning messages -- no longer printing '
         end if
         once = .true.
      else
         msg =   routine // ' warning '  // messag
      endif

      return

  end subroutine  errmsg

end module mo_errmsg
