!******************************************************************************
! DESCRIPTION:
!     Module to compute critical statistic for bootstrapping significance test.
!     SRS (simple random sampling) is used to choose a subset of individuals
!     (a sample) chosen from a larger set, which means each individual is chosen
!     randomly and entirely by chance, such that each individual has the same
!     probability of being chosen at any stage during the sampling process
!
! REVISION HISTORY:
!     Prototype 09/2016 by Shen Yu, SCC
!     Revised 09/2016 by Li Zhenkun, SCC
!
!******************************************************************************
  module critical_statistic_mod

     implicit none

     public       ::  upper_u0
     public       ::  sample
     private      ::  randnum
     private      ::  sort

     contains

!------------------------------------------------------------------------------
!     Public subroutine to find upper critical statistic u0 for significance test
!------------------------------------------------------------------------------
     subroutine upper_u0(n_resample, statistic, alpha, u0, ierr)

        implicit none

!------------------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(upper_u0)'
        integer, intent(in)                ::  n_resample
        real(4), intent(inout)             ::  statistic(n_resample)
        real(4), intent(inout)             ::  alpha
        real(4), intent(out)               ::  u0
        integer, intent(out)               ::  ierr

!------------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------------
        integer                            ::  q0 ! the number of upper extreme points in sample distribution which satisfies q0 / n_resample = alpha
        integer                            ::  n0 ! ascending-order-array statistic's subscript at which its sample distribution percentile is (1. - alpha)

        ierr = 0
        if ( alpha > 0.0 .and. alpha <= 0.25 ) then
           call sort(n_resample, statistic) ! sorting array statistic in ascending order
           q0 = int(alpha * real(n_resample)) ! finding the number of upper extreme points in sample distribution which satisfies q0 / n_resample = alpha
           n0 = n_resample - q0 + 1 ! finding ascending-order-array statistic's subscript at which its sample distributin percentile is (1. - alpha)
           u0 = statistic(n0)
           alpha = real(q0 + 1) / real(n_resample + 1) ! return revised more precise significant level alpha.
        else
           write( 6, '(2a)' ) subname, 'ERROR: alpha must be in (0, 0.25)'
           ierr = 1
           return
        end if

     end subroutine upper_u0

!------------------------------------------------------------------------------
!     Public subroutine to do SRS (simple random sample) from number1 to number2
!------------------------------------------------------------------------------
     subroutine sample(number1, number2, random_num)

        implicit none

!------------------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------------------
        integer, intent(in)                ::  number1
        integer, intent(in)                ::  number2
        integer, intent(out)               ::  random_num

!------------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------------
        integer                            ::  idum = -1

        random_num = min(number1, number2) + int(randnum(idum) * (abs(number2 - number1) + 1))

     end subroutine sample

!------------------------------------------------------------------------------
!     Private function to generate a uniform random deviate between 0.0 and 1.0
!     (exclusive of the endpoint values)
!------------------------------------------------------------------------------
     real function randnum(idum)

        implicit none

!------------------------------------------------------------------------------
!       dummy arguments
!------------------------------------------------------------------------------
        integer, intent(inout)             ::  idum
!      Call with idum a negative integer to initialize; thereafter, do not alter
!      idum between successive deviates in a sequence

!------------------------------------------------------------------------------
!       local arguments
!------------------------------------------------------------------------------
        integer, parameter                 ::  ia = 16807, im = 2147483647
        integer, parameter                 ::  iq = 127773, ir = 2836
        integer, parameter                 ::  ntab = 32, ndiv = 1 + (im - 1) / ntab
        real(4), parameter                 ::  am = 1. / im, eps = 1.2e-7, rnmx = 1. - eps
        integer                            ::  j, k, iv(ntab), iy
        save iv, iy
        data iv / ntab * 0 /, iy / 0 /

        if ( idum <= 0 .or. iy == 0 ) then
           idum = max( -idum, 1 )
           do j = ntab + 8, 1, -1
              k = idum / iq
              idum = ia * (idum - k * iq) - ir * k
              if ( idum < 0 ) idum = idum + im
              if ( j <= ntab ) iv(j) = idum
           end do
           iy = iv(1)
        end if

        k = idum / iq
        idum = ia * (idum - k * iq) - ir * k
        if ( idum < 0 ) idum = idum + im
        j = 1 + iy / ndiv
        iy = iv(j)
        iv(j) = idum
        randnum = min(am * iy, rnmx)

        return

     end function randnum

!------------------------------------------------------------------------------
!     Private subroutine to resort an array in ascending order
!------------------------------------------------------------------------------
     subroutine sort(n, x)

        implicit none

!------------------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------------------
        integer, intent(in)                ::  n
        real(4), intent(inout)             ::  x(n)

!------------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------------
        integer                            ::  i, j, l, ir
        real(4)                            ::  y

!------------------------------------------------------------------------------
        l = n / 2 + 1
        ir = n
        do
           if ( l > 1 ) then
              l = l - 1
              y = x(l)
           else
              y = x(ir)
              x(ir) = x(1)
              ir = ir - 1
              if ( ir == 1 ) then
                 x(1) = y
                 return
              end if
           end if
           i = l
           j = l + l
           do while ( j <= ir )
              if ( j < ir ) then
                 if ( x(j) < x(j+1) ) j = j + 1
              end if
              if ( y < x(j) ) then
                 x(i) = x(j)
                 i = j
                 j = j + j
              else
                 j = ir + 1
              end if
           end do
           x(i) = y
        end do

     end subroutine sort

  end module critical_statistic_mod
