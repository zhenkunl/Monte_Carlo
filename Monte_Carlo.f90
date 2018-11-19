!***********************************************************************
!  DESCRIPTION:
!       Main Program to verify dynamical analogue prediction using Monte
!       Carlo method
!
!  REVISION  HISTORY:
!       Prototype 09/2016 by Li Zhenkun, SCC
!       Revised 10/2018 by Li Zhenkun, SCC
!
!***********************************************************************
  program Monte_Carlo

     use params_mod
     use dynamical_analogue_mod
     use inout_mod
     use critical_statistic_mod
     use utils_mod

     implicit none

     integer, parameter                ::  stn_num = 523 ! observation station number
     real(4)                           ::  time_begin, time_end
     character(256)                    ::  progname, namelistfile
     character(256)                    ::  message
     integer                           ::  ierr
     integer                           ::  num_year
     integer                           ::  fct_year
     integer, allocatable              ::  samp_num(:)
     integer                           ::  samp_sgl_num
     integer                           ::  samp_max ! the maximum lenth of historical samples
     integer, allocatable              ::  samp_dates_valid(:, :, :)
     real(4), allocatable              ::  prec_fct_dap(:, :), prec_fct_dap_ens(:)
     real(4), allocatable              ::  prec_fct_random(:, :), prec_fct_random_ens(:)
     real(4), allocatable              ::  prec_obs(:, :)
     real(4), allocatable              ::  conc_dap(:, :), conc_random(:, :, :)
     real(4)                           ::  conc_dap_mean(2)
     real(4), allocatable              ::  conc_random_mean(:, :)
     real(4), allocatable              ::  acc_dap(:, :), acc_random(:, :, :)
     real(4)                           ::  acc_dap_mean(2)
     real(4), allocatable              ::  acc_random_mean(:, :)
     integer                           ::  iyear
     integer                           ::  isamp, inum
     integer                           ::  random_num
     real(4)                           ::  alpha0(4), u0(4)
     integer                           ::  iunit
     integer                           ::  status
     character(256)                    ::  filename
     character(2)                      ::  cmon, cday
     character(2)                      ::  lead_days_str, ave_days_str

     call cpu_time(time_begin)

!----------------------------------------------------------------------
!    Get the names of this program and namelist control file
!----------------------------------------------------------------------
     call getarg(0, progname)
     call getarg(1, namelistfile)
     call init_params(namelistfile, ierr)
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'Parameter initialization not completed'
        write( 6, '(a)' ) 'Usage : '
        write( 6, '(3a)' ) '        ', trim(progname), ' namelist'
        write( 6, '(a)' ) ' '
        write( 6, '(a)' ) 'Check argument and namelist syntax'
        stop
     end if

!----------------------------------------------------------------------
!    Check the length of leading days and average days are valid or not
!    according to the forecast period of NCEP CFSv2 (43 days), that is,
!    lead_days + ave_days should be NO more than 45 days
!----------------------------------------------------------------------
     if ( lead_days > 44 ) then
        write ( 6, '(a)' ) 'Leading days setting error, check namelist argument'
        stop
     end if
     if ( ave_days + lead_days > 45 ) then
        write ( 6, '(a)' ) 'Average days setting error, check namelist argument'
        stop
     end if

!----------------------------------------------------------------------
!    Check the domain settings are invalid or not
!----------------------------------------------------------------------
     if ( domain_west < 0 )   domain_west = domain_west + 360
     if ( domain_west >= 360 ) domain_west = domain_west - 360
     if ( domain_east < 0 )   domain_east = domain_east + 360
     if ( domain_east >= 360 ) domain_east = domain_east - 360
     if ( domain_west > domain_east .or. domain_south > domain_north ) then
        write( 6, '(a)' ) 'Domain settings error, check namelist argument'
        stop
     end if

     num_year = end_year - start_year + 1
     allocate( samp_num(num_year), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: samp_num'
        stop
     end if

     samp_sgl_num = samp_sgl_len / samp_interval
     samp_max = (samp_sgl_num * 2 + 1) * (end_year - start_year + 1)

     allocate( samp_dates_valid(3, samp_max, num_year), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: samp_dates_valid'
        stop
     end if

     allocate( prec_fct_dap(stn_num, anal_fct_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: prec_fct_dap'
        stop
     end if

     allocate( prec_fct_dap_ens(stn_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: prec_fct_dap_ens'
        stop
     end if

     allocate( prec_obs(stn_num, num_year), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: prec_obs'
        stop
     end if

     allocate( conc_dap(2, num_year), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: conc_dap'
        stop
     end if

     allocate( acc_dap(2, num_year), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: acc_dap'
        stop
     end if

     do fct_year = start_year, end_year

        iyear = fct_year - start_year + 1
        call get_prec_anom_direct(input_dir_prec, path_separator, stn_num, fct_year, fct_mon, fct_day, start_year, end_year, fct_mon, fct_day, lead_days, ave_days, prec_obs(:, iyear), ierr)
        if ( ierr /= 0 ) then
           write( 6, '(a)' ) 'ERROR: can NOT compute precipitation anomaly percentage'
           stop
        end if

        write ( message, '(a, i4, 1x, i2, 1x, i2)' ) 'Begin to compute dynamical analogue prediction for ', fct_year, fct_mon, fct_day
        call print_message(100, debug_level, message)
        call dynamical_analogue(fct_year, fct_mon, fct_day, lead_days, ave_days, domain_west, domain_east, domain_south, domain_north, start_year, end_year, samp_interval, &
                                samp_sgl_num, anal_fct_num, var_num, var_names, input_dir_curt, input_dir_hist, input_dir_prec, path_separator, stn_num, samp_max, samp_num(iyear), &
                                samp_dates_valid(:, :, iyear), prec_fct_dap, debug_level, ierr)

        prec_fct_dap_ens = sum(prec_fct_dap, dim=2) / real(anal_fct_num)


        call get_prec_conc(stn_num, prec_fct_dap(:,1), prec_obs(:, iyear), conc_dap(1, iyear))
        call get_prec_conc(stn_num, prec_fct_dap_ens, prec_obs(:, iyear), conc_dap(2, iyear))

        call get_prec_acc(stn_num, prec_fct_dap(:,1), prec_obs(:, iyear), acc_dap(1, iyear))
        call get_prec_acc(stn_num, prec_fct_dap_ens, prec_obs(:, iyear), acc_dap(2, iyear))

     end do

     conc_dap_mean = sum(conc_dap, dim=2) / real(num_year)
     acc_dap_mean = sum(acc_dap, dim=2) / real(num_year)

     deallocate( prec_fct_dap )
     deallocate( prec_fct_dap_ens )
     deallocate( conc_dap )
     deallocate( acc_dap )

     allocate( prec_fct_random(stn_num, anal_fct_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: prec_fct_random'
        stop
     end if

     allocate( prec_fct_random_ens(stn_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: prec_fct_random_ens'
        stop
     end if

     allocate( conc_random(2, num_year, resample_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: conc_random'
        stop
     end if

     allocate( conc_random_mean(2, resample_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: conc_random_mean'
        stop
     end if

     allocate( acc_random(2, num_year, resample_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: acc_random'
        stop
     end if

     allocate( acc_random_mean(2, resample_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: acc_random_mean'
        stop
     end if

     do inum = 1, resample_num
        write( message, '(a, i, a)' ) 'Begin the ', inum, 'th simple random sampling'
        call print_message(100, debug_level, message)

        do fct_year = start_year, end_year
           write( message, '(a, i4)' ) '   Begin to choose random similar sample prediction for ', fct_year
           call print_message(200, debug_level, message)
           iyear = fct_year - start_year + 1
           do isamp = 1, anal_fct_num
              call sample(1, samp_num(iyear), random_num)
              call get_prec_anom_direct(input_dir_prec, path_separator, stn_num, samp_dates_valid(1, random_num, iyear), samp_dates_valid(2, random_num, iyear), &
                                        samp_dates_valid(3, random_num, iyear), start_year, end_year, fct_mon, fct_day, lead_days, ave_days, prec_fct_random(:, isamp), ierr)
              if ( ierr /= 0 ) then
                 write( 6, '(a)' ) 'ERROR: can NOT compute precipitation anomaly percentage'
                 stop
              end if
           end do

           prec_fct_random_ens = sum(prec_fct_random, dim=2)/ real(anal_fct_num)

           call get_prec_conc(stn_num, prec_fct_random(:,1), prec_obs(:, iyear), conc_random(1, iyear, inum))
           call get_prec_conc(stn_num, prec_fct_random_ens, prec_obs(:, iyear), conc_random(2, iyear, inum))

           call get_prec_acc(stn_num, prec_fct_random(:,1), prec_obs(:, iyear), acc_random(1, iyear, inum))
           call get_prec_acc(stn_num, prec_fct_random_ens, prec_obs(:, iyear), acc_random(2, iyear, inum))

        end do

        conc_random_mean = sum(conc_random, dim=2) / real(num_year)
        acc_random_mean = sum(acc_random, dim=2) / real(num_year)

     end do

     alpha0 = alpha
     call upper_u0(resample_num, conc_random_mean(1, :), alpha0(1), u0(1), ierr)
     if ( ierr /= 0 ) then
       write( 6, '(a)' ) 'ERROR: can NOT find upper critical statistic u0 for significance test'
       stop
     end if
     call upper_u0(resample_num, conc_random_mean(2, :), alpha0(2), u0(2), ierr)
     if ( ierr /= 0 ) then
       write( 6, '(a)' ) 'ERROR: can NOT find upper critical statistic u0 for significance test'
       stop
     end if
     call upper_u0(resample_num, acc_random_mean(1, :), alpha0(3), u0(3), ierr)
     if ( ierr /= 0 ) then
       write( 6, '(a)' ) 'ERROR: can NOT find upper critical statistic u0 for significance test'
       stop
     end if
     call upper_u0(resample_num, acc_random_mean(2, :), alpha0(4), u0(4), ierr)
     if ( ierr /= 0 ) then
       write( 6, '(a)' ) 'ERROR: can NOT find upper critical statistic u0 for significance test'
       stop
     end if

     iunit = 100
     cmon = char(fct_mon/10 + 48)//char(mod(fct_mon, 10) + 48)
     cday = char(fct_day/10 + 48)//char(mod(fct_day, 10) + 48)
     lead_days_str = char(lead_days/10 + 48)//char(mod(lead_days, 10) + 48)
     ave_days_str = char(ave_days/10 + 48)//char(mod(ave_days, 10) + 48)
     filename = trim(output_dir)//path_separator//cmon//cday//'-lead-'//lead_days_str//'-for-'//ave_days_str//'-forecast-MonteCarlo.txt'
     open( unit = iunit, file = trim(filename), form = 'formatted', iostat = status )
     if ( status /= 0 ) then
        write( 6, '(2a)' ) 'ERROR: can NOT rightly open file ', trim(filename)
        stop
     end if

     write( iunit, '(<resample_num>(1x, f8.5))' ) conc_random_mean(1, :)
     write( iunit, '(3(1x,f8.5))' ) alpha0(1), u0(1), conc_dap_mean(1)
     write( iunit, '(<resample_num>(1x, f8.5))' ) conc_random_mean(2, :)
     write( iunit, '(3(1x,f8.5))' ) alpha0(2), u0(2), conc_dap_mean(2)

     write( iunit, '(<resample_num>(1x, f8.5))' ) acc_random_mean(1, :)
     write( iunit, '(3(1x,f8.5))' ) alpha0(3), u0(3), acc_dap_mean(1)
     write( iunit, '(<resample_num>(1x, f8.5))' ) acc_random_mean(2, :)
     write( iunit, '(3(1x,f8.5))' ) alpha0(4), u0(4), acc_dap_mean(2)

     close( iunit )

     deallocate( samp_num )
     deallocate( samp_dates_valid )
     deallocate( prec_fct_random )
     deallocate( prec_fct_random_ens )
     deallocate( prec_obs )
     deallocate( conc_random )
     deallocate( conc_random_mean )
     deallocate( acc_random )
     deallocate( acc_random_mean )

     write( 6, '(a)' ) '=============== SUCCESSFUL TERMINATION OF MONTE CARLO TEST ==============='
     call cpu_time(time_end)
     write( 6, '(a, f9.3, a)' ) '=============== compute time (seconds)      =   ', time_end - time_begin, ' ==============='

  end program Monte_Carlo
