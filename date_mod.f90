!***********************************************************************
!  DESCRIPTION:
!       Module to process date related issues
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!       Revised 10/2018 by Li Zhenkun, SCC
!
!***********************************************************************
  module date_mod

     implicit none

     public     ::  new_date
     private    ::  nfeb

     contains

!----------------------------------------------------------------------
!    Public subroutine to calculate the new date for a given start
!    date and interval
!----------------------------------------------------------------------
     subroutine new_date(iyr, imo, idy, incre, jyr, jmo, jdy)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer, intent(in)         ::  iyr, imo, idy, incre
        integer, intent(out)        ::  jyr, jmo, jdy

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        integer                     ::  days_month(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        integer                     ::  nday, i

        jyr  = iyr
        jmo  = imo
        jdy  = idy

        if ( incre == 0) return

        nday = 0
        do i = 1, imo
           if ( i == 2 ) days_month(i) = nfeb(iyr)
           nday = nday + days_month(i)
        end do
        nday = nday - (days_month(imo) - idy)

        nday = nday + incre

        if ( incre > 0 ) then
           do while ( nday > 337+nfeb(jyr) )
              nday = nday - ( 337 + nfeb(jyr) )
              jyr = jyr + 1
           end do
        else
           do while ( nday <= 0 )
              jyr = jyr - 1
              nday = nday + 337 + nfeb(jyr)
           end do
        end if

        days_month(2) = nfeb(jyr)
        jmo = 1
        do while ( nday > days_month(jmo) )
           nday = nday - days_month(jmo)
           jmo = jmo + 1
        end do
        jdy = nday

     end subroutine new_date

!----------------------------------------------------------------------
!    Private function to calculate the number of day in February for a
!    given year
!----------------------------------------------------------------------
     integer function nfeb(iyear)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer           ::  iyear

        if ( (mod(iyear,4) == 0 .and. mod(iyear,100) /= 0) .or. mod(iyear,400) == 0 ) then
           nfeb = 29
        else
           nfeb = 28
        end if

        return

     end function nfeb

  end module date_mod
