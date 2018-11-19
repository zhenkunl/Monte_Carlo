!***********************************************************************
!  DESCRIPTION:
!       Module to predict short-term climate using dynamical analogue
!       method based on NCEP CFSv2 operational products
!
!  PRECONDITIONS REQUIRED:
!       NCEP CFSv2 raw data should be processed on Unix-like platform
!       or in Windows/Cygwin environment first
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!       Revised 09/2016 by Li Zhenkun, SCC
!       Revised 10/2018 by Li Zhenkun, SCC
!
!***********************************************************************
  module dynamical_analogue_mod

     implicit none

     public       ::  dynamical_analogue

     contains

!------------------------------------------------------------------
!    Public subroutine to predict short-term climate using dynamical
!    analogue method
!------------------------------------------------------------------
     subroutine dynamical_analogue(fct_year, fct_mon, fct_day, lead_days, ave_days, domain_west, domain_east, domain_south, domain_north, start_year, end_year, samp_interval, &
                                   samp_sgl_num, anal_fct_num, var_num, var_names, input_dir_curt, input_dir_hist, input_dir_prec, path_separator, stn_num, samp_max, samp_num, &
                                   samp_dates_valid, prec_fct, debug_level, ierr)

        use inout_mod
        use eof_mod
        use date_mod, only :  new_date
        use utils_mod

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer, intent(in)               ::  fct_year, fct_mon, fct_day
        integer, intent(in)               ::  lead_days, ave_days
        integer, intent(in)               ::  domain_west, domain_east, domain_south, domain_north
        integer, intent(in)               ::  start_year, end_year
        integer, intent(in)               ::  samp_interval, samp_sgl_num
        integer, intent(inout)            ::  anal_fct_num ! selected forecast number
        integer, intent(in)               ::  var_num ! for further research purpose
        character(len=*), intent(in)      ::  var_names(var_num)
        character(len=*), intent(in)      ::  input_dir_curt, input_dir_hist, input_dir_prec, path_separator
        integer, intent(in)               ::  stn_num
        integer, intent(in)               ::  samp_max
        integer, intent(out)              ::  samp_num
        integer, intent(out)              ::  samp_dates_valid(3, samp_max)
        real(4), intent(out)              ::  prec_fct(stn_num, anal_fct_num)
        integer, intent(in)               ::  debug_level
        integer, intent(out)              ::  ierr ! error identifier

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        character(len=*), parameter       ::  subname = '(dynamical_analogue)'
        integer                           ::  nlon, nlat
        integer                           ::  grid_num
        real(4), allocatable              ::  curt(:), hist_single(:), hist_valid(:, :), hist(:, :)
        character(256)                    ::  message
        integer                           ::  incre
        integer                           ::  year, mon, day
        integer, allocatable              ::  samp_dates(:, :)
        integer                           ::  days_month(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        real(4), allocatable              ::  mean(:), sdev(:)
        integer                           ::  abnorm_num
        real(4)                           ::  minima = 1.e-25
        real(4)                           ::  undef = -999.0
        integer                           ::  min_dim
        integer                           ::  trans_opt = 2 !trans_opt = 0 for raw, 1 for anomaly, 2 for normalized anomoly
        real(4), allocatable              ::  eigenvalue(:)
        real(4), allocatable              ::  eigenvector(:, :)
        real(4), allocatable              ::  pc(:, :)
        real(4)                           ::  eigenvalue_sum
        integer                           ::  lead_num
        real(4), allocatable              ::  proj_coeff_hist(:, :)
        real(4), allocatable              ::  proj_coeff_curt(:)
        real(4), allocatable              ::  similar_val(:)
        integer                           ::  iyear
        integer                           ::  igrid, isamp, inum

        nlon = domain_east - domain_west + 1
        nlat = domain_north - domain_south + 1
        grid_num = nlon * nlat * var_num

!----------------------------------------------------------------------
!    Get varibles of current forecast from pre-processed binary file
!----------------------------------------------------------------------
        message = subname//'Step 1: acquire current forecast data ...'
        call print_message(200, debug_level, message)

        allocate( curt(grid_num), stat = ierr )
        if ( ierr /= 0 ) then
          write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: curt'
          return
        end if

        call get_cfs_field(input_dir_curt, path_separator, var_num, var_names(1:var_num), fct_year, fct_mon, fct_day, &
                           lead_days, ave_days, domain_west, domain_south, nlon, nlat, curt, ierr)
        if ( ierr /= 0 ) then
           write( 6, '(2a, i4, 1x, i2, 1x, i2, a)' ) subname, 'ERROR: current forecast data of ', fct_year, fct_mon, fct_day, ' is incomplete'
           return
        end if

!----------------------------------------------------------------------
!    Get varibles of historical forecast from pre-processed binary file
!----------------------------------------------------------------------
        message = subname//'Step 2: acquire historical forecast data ...'
        call print_message(200, debug_level, message)

        allocate( hist_single(grid_num), stat = ierr )
        if ( ierr /= 0 ) then
          write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: hist_single'
          return
        end if

        allocate( hist_valid(grid_num, samp_max), stat = ierr )
        if ( ierr /= 0 ) then
          write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: hist_valid'
          return
        end if

        samp_num = 0
        do iyear = start_year, end_year
           if ( iyear == fct_year ) cycle
           do isamp = -samp_sgl_num, samp_sgl_num

              incre = isamp * samp_interval
              call new_date(iyear, fct_mon, fct_day, incre, year, mon, day)
              if ( year < start_year .or. year > end_year ) cycle
              call get_cfs_field(input_dir_hist, path_separator, var_num, var_names(1:var_num), year, mon, day, &
                                 lead_days, ave_days, domain_west, domain_south, nlon, nlat, hist_single, ierr)
              if ( ierr /= 0 ) then
                 write( message, '(a, i4, 1x, i2, 1x, i2, a)' ) 'WARNING: forecast data of ', year, mon, day, ' is incomplete'
                 call print_message(200, debug_level, message)
                 cycle
              else
                 samp_num = samp_num + 1
                 hist_valid(:, samp_num) = hist_single
                 samp_dates_valid(1, samp_num) = year
                 samp_dates_valid(2, samp_num) = mon
                 samp_dates_valid(3, samp_num) = day
              end if

           end do
        end do
        write( message, '(a, i)' ) 'Valid number of forecast samples is ', samp_num
        call print_message(200, debug_level, message)

        if ( samp_num < anal_fct_num ) anal_fct_num = samp_num

        deallocate( hist_single )

        allocate( hist(grid_num, samp_num), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: hist'
           return
        end if

        allocate( samp_dates(3, samp_num), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: samp_dates'
           return
        end if

        hist(:, :) = hist_valid(:, 1:samp_num)
        samp_dates(:, :) = samp_dates_valid(:, 1:samp_num)

        deallocate( hist_valid )

!----------------------------------------------------------------------
!    Normalize current forecast data using mean and standard deviation
!    of historical forecast samples
!----------------------------------------------------------------------
        message = subname//'Step 3: normalize current forecast data ...'
        call print_message(200, debug_level, message)

        allocate( mean(grid_num), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: mean'
           return
        end if

        allocate( sdev(grid_num), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: sdev'
           return
        end if

        do igrid = 1, grid_num
           call statistic(samp_num, hist(igrid, :), mean(igrid), sdev(igrid))
        end do

        abnorm_num = count(sdev <= minima)
        if ( abnorm_num /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: Undefind value found in matrix standrad deviation'
           ierr = 1
           return
        end if
        curt = (curt - mean) / sdev

        deallocate( mean, sdev )

!----------------------------------------------------------------------
!    Do Multivariate EOF decomposition to historical forecast to obtain
!    eigenvalues and eigenvectors
!----------------------------------------------------------------------
        message = subname//'Step 4: do EOF decomposition ...'
        call print_message(200, debug_level, message)

        min_dim = min( grid_num, samp_num )
        allocate( eigenvalue(min_dim), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: eigenvalue'
           return
        end if

        allocate( eigenvector(grid_num, min_dim), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: eigenvector'
           return
        end if

        allocate( pc(min_dim, samp_num), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: pc'
           return
        end if

        call eof(grid_num, samp_num, min_dim, trans_opt, undef, hist, eigenvalue, eigenvector, pc, ierr)

        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: EOF decomposition can NOT be rightly done'
           return
        end if

        deallocate( hist )

!----------------------------------------------------------------------
!    Seek the number of the first several leading EOFs whose cumulative
!    variance contribute rate excceed a certain level
!----------------------------------------------------------------------
        message = subname//'Step 5: seek suitable number of EOFs ...'
        call print_message(200, debug_level, message)

        eigenvalue_sum = sum(eigenvalue)
        lead_num = 10
        write( message, '(a, i2, a)' ) 'The cumulative variance contribute rates of the first ', lead_num, ' leading EOFs are: '
        call print_message(200, debug_level, message)
        do inum = 1, lead_num
           write( message, '(a, i2, 2x, f6.3, a)' ) '        ', inum, 100.0 * sum(eigenvalue(1:inum)) / eigenvalue_sum, '%'
           call print_message(200, debug_level, message)
        end do

        deallocate( eigenvalue )

        allocate( proj_coeff_hist(lead_num, samp_num), stat = ierr )
        if ( ierr /= 0 ) then
           write ( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: proj_coeff_hist'
           return
        end if

        allocate( proj_coeff_curt(lead_num), stat = ierr )
        if ( ierr /= 0 ) then
           write ( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: proj_coeff_curt'
           return
        end if

        proj_coeff_hist(:, :) = pc(1:lead_num, :)

        deallocate( pc )

!----------------------------------------------------------------------
!    Calculate projection coefficients using dot product of current
!    forecast and selected eigenvectors after current forecast is
!    normalized
!----------------------------------------------------------------------
        message = subname//'Step 6: calculate projection coefficients ...'
        call print_message(200, debug_level, message)

        do inum = 1, lead_num
           proj_coeff_curt(inum) = dot_product(eigenvector(:, inum), curt(:))
        end do

        deallocate( eigenvector )
        deallocate( curt )

!----------------------------------------------------------------------
!    Calculate similarity coefficients between current forecast and
!    historical forecasts
!----------------------------------------------------------------------
        message = subname//'Step 7: calculate similarity coefficients ...'
        call print_message(200, debug_level, message)

        allocate( similar_val(samp_num), stat = ierr )
        if ( ierr /= 0 ) then
           write( 6, '(2a)' ) subname, 'ERROR: can NOT allocate memory space to varible: similar_val'
           return
        end if

        do isamp = 1, samp_num
           call simi_disc_val(lead_num, proj_coeff_curt, proj_coeff_hist(:, isamp), similar_val(isamp))
        end do

        deallocate( proj_coeff_curt )
        deallocate( proj_coeff_hist )

!----------------------------------------------------------------------
!    Reorder similarity coefficients in ascending order to obtain specific
!    dates of the most similar historical forecasts
!----------------------------------------------------------------------
        message = subname//'Step 8: reorder similarity coefficients ...'
        call print_message(200, debug_level, message)

        call sort(samp_num, 3, similar_val, samp_dates)
        write( message, '(a, i3, a)' ) 'The first ', anal_fct_num, ' most similar discrete values are: '
        call print_message(200, debug_level, message)
        do isamp = 1, anal_fct_num
           write( message, '(a, i2, 1x, f)' ) '        ', isamp, similar_val(isamp)
           call print_message(200, debug_level, message)
        end do

        write( message, '(a, i3, a)' ) 'The first ', anal_fct_num, ' most similar forecasts dates are: '
        call print_message(200, debug_level, message)
        do isamp = 1, anal_fct_num
           write( message, '(a, i2, 2x, i4, 1x, i2, 1x, i2)' ) '        ', isamp, samp_dates(1, isamp), samp_dates(2, isamp), samp_dates(3, isamp)
           call print_message(200, debug_level, message)
        end do

        deallocate( similar_val )

!----------------------------------------------------------------------
!    Calculate precipitation anomaly percentage of the most similar
!    historical forecasts
!----------------------------------------------------------------------
        message = subname//'Step 9: calculate precipitation anomaly percentage ...'
        call print_message(200, debug_level, message)

        do isamp = 1, anal_fct_num
           call get_prec_anom_direct(input_dir_prec, path_separator, stn_num, samp_dates(1, isamp), samp_dates(2, isamp), &
                                     samp_dates(3, isamp), start_year, end_year, fct_mon, fct_day, lead_days, ave_days, prec_fct(:, isamp), ierr)
           if ( ierr /= 0 ) then
              write( 6, '(2a)' ) subname, 'ERROR: can NOT compute precipitation anomaly percentage'
              return
           end if
        end do

        deallocate( samp_dates )

     end subroutine dynamical_analogue

  end module dynamical_analogue_mod
