!***********************************************************************
!  DESCRIPTION:
!       Module for some useful utilities
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!
!***********************************************************************
  module utils_mod

     implicit none

     public      ::  simi_disc_val
     public      ::  sort
     public      ::  get_prec_conc
     public      ::  get_prec_acc
     public      ::  cor_coeff
     public      ::  print_message

     contains

!------------------------------------------------------------------
!     Public subroutine to compute Similar Discrete Values between
!     two vectors
!------------------------------------------------------------------
     subroutine simi_disc_val(n, x, y, value)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n
        real(4), intent(in)                ::  x(n), y(n)
        real(4), intent(out)               ::  value

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        integer                            ::  i
        real(4)                            ::  dx(n)
        real(4)                            ::  s, d, e

        do i = 1, n
           dx(i) = x(i) - y(i)
        end do

        e = 0.
        d = 0.
        do i = 1, n
           e = e + dx(i)
           d = d + abs(dx(i))
        end do
        e = e / real(n)
        d = d / real(n)

        s = 0.
        do i = 1, n
           s = s + abs(dx(i) - e)
        end do
        s = s / real(n)

        value = (s + d) / 2.0

     end subroutine simi_disc_val

!------------------------------------------------------------------
!     Public subroutine to resort an array in ascending order and
!     its corresponding two-dimension array
!------------------------------------------------------------------
     subroutine sort(n, m, x, y)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n, m
        real(4), intent(inout)             ::  x(n)
        integer, intent(inout)             ::  y(m, n)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        integer                            ::  l, ir
        real(4)                            ::  rra, rrb(3)
        integer                            ::  i, j

        l = n/2 + 1
        ir = n
        do
           if ( l > 1 ) then
              l = l - 1
              rra = x(l)
              rrb(:) = y(:, l)
           else
              rra = x(ir)
              rrb(:) = y(:, ir)
              x(ir) = x(1)
              y(:, ir) = y(:, 1)
              ir = ir - 1
              if ( ir == 1 ) then
                 x(1) = rra
                 y(:, 1) = rrb(:)
                 return
              end if
           end if
           i = l
           j = l + l
           do while ( j <= ir )
              if ( j < ir ) then
                 if ( x(j) < x(j+1) ) j = j + 1
              end if
              if ( rra < x(j) ) then
                 x(i) = x(j)
                 y(:, i) = y(:, j)
                 i = j
                 j = j + j
              else
                 j = ir + 1
              endif
           end do
           x(i) = rra
           y(:, i) = rrb(:)
        end do

     end subroutine sort

!------------------------------------------------------------------
!     Public subroutine to calculate precipitation concordance
!------------------------------------------------------------------
     subroutine get_prec_conc(n, x, y, conc)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n
        real(4), intent(in)                ::  x(n), y(n)
        real(4), intent(out)               ::  conc

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        integer                            ::  i, num

        num = 0
        do i = 1, n
           if ( x(i) * y(i) > 0. ) num = num + 1
        end do

        conc = num / real(n)

     end subroutine get_prec_conc

!------------------------------------------------------------------
!     Public subroutine to calculate precipitation ACC
!------------------------------------------------------------------
     subroutine get_prec_acc(n, x, y, acc)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n
        real(4), intent(in)                ::  x(n), y(n)
        real(4), intent(out)               ::  acc

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        integer                            ::  i
        real(4)                            ::  ave_var1, ave_var2
        real(4)                            ::  dev_var1, dev_var2
        real(4)                            ::  denominator, numerator

        ave_var1 = 0.
        ave_var2 = 0.
        do i = 1, n
           ave_var1 = ave_var1 + x(i)
           ave_var2 = ave_var2 + y(i)
        end do
        ave_var1 = ave_var1 / real(n)
        ave_var2 = ave_var2 / real(n)

        dev_var1 = 0.
        dev_var2 = 0.
        do i = 1, n
           dev_var1 = dev_var1 + ( x(i) - ave_var1 ) ** 2
           dev_var2 = dev_var2 + ( y(i) - ave_var2 ) ** 2
        end do
        dev_var1 = sqrt(dev_var1)
        dev_var2 = sqrt(dev_var2)

        numerator = 0.
        do i = 1, n
           numerator = numerator + (x(i) - ave_var1) * (y(i) - ave_var2)
        end do

        denominator = dev_var1 * dev_var2

        acc = numerator / denominator

     end subroutine get_prec_acc

!------------------------------------------------------------------
!     Public subroutine to calculate correlation coefficient between
!     two three-dimension meteorological fields
!------------------------------------------------------------------

     subroutine cor_coeff(nlon, nlat, nsamp, var1, var2, cor)

        implicit none
!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  nlon, nlat, nsamp
        real(4), intent(in)                ::  var1(nlon, nlat, nsamp)
        real(4), intent(in)                ::  var2(nlon, nlat, nsamp)
        real(4), intent(out)               ::  cor(nlon, nlat)
!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        integer                            ::  ilon, ilat, isamp
        real(4)                            ::  ave_var1, ave_var2
        real(4)                            ::  dev_var1, dev_var2
        real(4)                            ::  denominator, numerator

        do ilon = 1, nlon
           do ilat = 1, nlat
              ave_var1 = 0.
              ave_var2 = 0.
              do isamp = 1, nsamp
                 ave_var1 = ave_var1 + var1(ilon, ilat, isamp)
                 ave_var2 = ave_var2 + var2(ilon, ilat, isamp)
              end do
              ave_var1 = ave_var1/nsamp
              ave_var2 = ave_var2/nsamp

              dev_var1 = 0.
              dev_var2 = 0.
              do isamp = 1, nsamp
                 dev_var1 = dev_var1 + (var1(ilon, ilat, isamp) - ave_var1)**2
                 dev_var2 = dev_var2 + (var2(ilon, ilat, isamp) - ave_var2)**2
              end do
              dev_var1 = sqrt(dev_var1)
              dev_var2 = sqrt(dev_var2)

              numerator = 0.
              do isamp = 1, nsamp
                 numerator = numerator + (var1(ilon, ilat, isamp) - ave_var1)*(var2(ilon, ilat, isamp) - ave_var2)
              end do

              denominator = dev_var1 * dev_var2

              cor(ilon, ilat) = numerator/denominator

           end do
        end do

     end subroutine cor_coeff

!------------------------------------------------------------------
!     Public subroutine to print debug message
!------------------------------------------------------------------
     subroutine print_message(level, debug_level, message)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  level, debug_level
        character(len=*), intent(in)       ::  message

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------

        if ( level <= debug_level ) then
           write( 6, '(a)' ) message
        end if

     end subroutine print_message

  end module utils_mod
