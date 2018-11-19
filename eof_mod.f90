!***********************************************************************
!  DESCRIPTION:
!       Module applies the EOF approach to analysis meteorological field
!       of two dimensions
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Shen Yu, SCC
!       Revised 12/2015 by Li Zhenkun, SCC
!
!***********************************************************************
  module eof_mod

     implicit none

     public      ::  statistic
     public      ::  eof

     contains

!----------------------------------------------------------------------
!    Public subroutine to compute mean and standard deviation statistic
!    of a one-dimension array
!----------------------------------------------------------------------
     subroutine statistic(n, array, mean, sdev)

        implicit none
!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n
        real(4), intent(in)                ::  array(n)
        real(4), intent(out)               ::  mean, sdev

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        integer                            ::  i

        mean = 0.
        do i = 1, n
           mean = mean + array(i)
        end do
        mean = mean / real(n)

        sdev = 0.
        do i = 1, n
            sdev = sdev + ( array(i) - mean ) ** 2
        end do
        sdev = sqrt(sdev / real(n))

     end subroutine statistic

!----------------------------------------------------------------------
!    Public subroutine to do EOF decomposition
!----------------------------------------------------------------------
     subroutine eof(grid_num, samp_num, min_dim, trans_opt, undef, field, eigenvalue, eigenvector, pc, ierr)

        use eigen_mod

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  grid_num, samp_num, min_dim
        integer, intent(in)                ::  trans_opt !trans_opt = 0 for raw, 1 for anomaly, 2 for normalized anomoly
        real(4), intent(in)                ::  undef
        real(4), intent(inout)             ::  field(grid_num, samp_num)
        real(4), intent(out)               ::  eigenvalue(min_dim)
        real(4), intent(out)               ::  eigenvector(grid_num, min_dim)
        real(4), intent(out)               ::  pc(min_dim, samp_num)
        integer, intent(out)               ::  ierr

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(eof)'
        integer, parameter                 ::  opt = 3 !opt = 1, 2 find eigenvalues only, opt = 3, 4 find both eigenvalues and its eigenvectors.
        real(4)                            ::  symm_matrix(min_dim, min_dim), symm_eigenvector(min_dim, min_dim)
        real(4)                            ::  mean(grid_num), sdev(grid_num)
        real(4)                            ::  lamda(min_dim)
        integer                            ::  abnorm_num
        real(4)                            ::  minima = 1.e-25
        integer                            ::  i, j

!----------------------------------------------------------------------
!    Check whether there are undef values in matrix field
!----------------------------------------------------------------------
        ierr = 0
        abnorm_num = count(field == undef)
        if ( abnorm_num /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: Undefind value found in matrix'
           ierr = 1
           return
        end if

!----------------------------------------------------------------------
!    Compute mean and standard deviation of matrix field for second
!    dimension samp_num
!----------------------------------------------------------------------
        do i = 1, grid_num
            call statistic(samp_num, field(i, :), mean(i), sdev(i))
        end do

!----------------------------------------------------------------------
!    Matrix transform according to the given parameter trans_opt
!----------------------------------------------------------------------
        if ( trans_opt == 1 ) then  ! anomaly transform

           forall ( i = 1:samp_num ) field(1:grid_num, i) = field(1:grid_num, i) - mean(1:grid_num)

        else if ( trans_opt == 2 ) then  ! standardized transform

           abnorm_num = count(sdev <= minima)  ! checking whether there are zeros in sdev
           if ( abnorm_num /= 0 ) then
              write( 6, '(2a)' ) subname, 'ERROR: Values very close to Zero found in matrix standrad deviation'
              ierr = 1
              return
           end if
           forall ( i = 1:samp_num ) field(1:grid_num, i) = ( field(1:grid_num, i) - mean(1:grid_num) ) / sdev(1:grid_num)

        end if

!----------------------------------------------------------------------
!    Transpose matrix if grid_num > samp_num (i.e. space dimension >
!    time dimension
!----------------------------------------------------------------------
        symm_matrix = 0.0
        if ( grid_num > samp_num ) then

           symm_matrix = matmul(transpose(field), field)
           symm_matrix = symm_matrix / grid_num

        else

           symm_matrix = matmul(field, transpose(field))
           symm_matrix = symm_matrix / samp_num

        end if

!----------------------------------------------------------------------
!    Call subroutine rs() to obtain all eigenvalues & corresponding
!    eigenvectors of symmetric matrix
!----------------------------------------------------------------------
        call rs(min_dim, min_dim, symm_matrix, eigenvalue, opt, symm_eigenvector, ierr)
        if ( ierr /= 0 ) then
           write( 6, '(2a, i, a)' ) subname, 'ERROR: can NOT evaluate the ', ierr, ' eigenvalue'
           return
        end if

!----------------------------------------------------------------------
!    Reorder the eigenvalues in a descending order as well as corres-
!    ponding eigenvectors
!----------------------------------------------------------------------
        eigenvalue(1:min_dim) = eigenvalue(min_dim:1:-1)
        symm_eigenvector(:, 1:min_dim) = symm_eigenvector(:, min_dim:1:-1)

        if ( grid_num > samp_num ) then ! if grid_num > samp_num find eigenvalue and its eigenvector

           eigenvalue = (grid_num / samp_num) * eigenvalue  ! eigenvalue
           eigenvector = matmul(field, symm_eigenvector) ! eigenvector
           forall( i = 1:min_dim ) lamda(i) =  sqrt(dot_product(eigenvector(:, i), eigenvector(:, i)))
           forall(i = 1:min_dim) eigenvector(:, i) = eigenvector(:, i) / lamda(i) ! normalized treatment.

        else

           eigenvector = symm_eigenvector

        end if

!----------------------------------------------------------------------
!    Obtain PCs of matrix field
!----------------------------------------------------------------------
        pc = matmul(transpose(eigenvector), field)

     end subroutine eof

  end module eof_mod
