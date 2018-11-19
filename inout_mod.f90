!***********************************************************************
!  DESCRIPTION:
!       Module for file I/O
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!       Revised 10/2018 by Li Zhenkun, SCC
!
!***********************************************************************
  module inout_mod

     implicit none

     public      ::  get_cfs_field
     public      ::  get_prec_anom_direct
     public      ::  get_prec_anom
     private     ::  get_prec

     contains

!----------------------------------------------------------------------
!     Public subroutine to get CFSv2 forecast in a given date, leading
!     days and average days
!----------------------------------------------------------------------
     subroutine get_cfs_field(input_dir, path_separator, var_num, var_names, year, mon, day, lead_days, ave_days, west, south, nlon, nlat, var, ierr)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        character(len=*), intent(in)       ::  input_dir
        character(len=*), intent(in)       ::  path_separator
        integer, intent(in)                ::  var_num
        character(len=*), intent(in)       ::  var_names(var_num)
        integer, intent(in)                ::  year, mon, day
        integer, intent(in)                ::  lead_days, ave_days
        integer, intent(in)                ::  west, south
        integer, intent(in)                ::  nlon, nlat
        real(4), intent(out)               ::  var(nlon*nlat*var_num)
        integer, intent(out)               ::  ierr

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(get_cfs_field)'
        character(4)                       ::  cyear
        character(2)                       ::  cmon, cday
        integer                            ::  istart, iend, jstart, jend
        character(256)                     ::  filename
        integer                            ::  iunit
        logical                            ::  lexist
        integer                            ::  status
        real(4)                            ::  var_whole(360, 181, ave_days)
        real(4)                            ::  var_region(nlon, nlat, ave_days)
        real(4)                            ::  var_local(nlon, nlat, var_num)
        integer                            ::  ivar, iday, irec

        write (cyear, '(i4)') year
        cmon = char(mon/10 + 48)//char(mod(mon, 10) + 48)
        cday = char(day/10 + 48)//char(mod(day, 10) + 48)

        istart = west + 1
        iend = istart + nlon -1
        jstart = south + 90 + 1
        jend = jstart + nlat - 1

        ierr = 0
        iunit = 101
        do ivar = 1, var_num

           filename = trim(input_dir)//path_separator//trim(var_names(ivar))//'/'//trim(var_names(ivar))//'-'//cyear//cmon//cday//'-cfs-forecast.dat'

           inquire( file = filename, exist = lexist )
           if ( .not. lexist ) then
              write( 6, '(4a)' ) subname, 'WARNING: ', trim(filename) , ' is NOT available'
              ierr = 1
              return
           end if

           open( unit = iunit, file = filename, form = 'unformatted', access = 'direct', recl = 360*181, iostat = status )
           if ( status /= 0 ) then
              write( 6, '(3a)' ) subname, 'WARNING: can NOT rightly open file ', trim(filename)
              ierr = 1
              return
           end if

           irec = lead_days - 1
           do iday = 1, ave_days
              read(iunit, rec = irec) var_whole(:, :, iday)
              irec = irec + 1
           end do
           var_region(:, :, :) = var_whole(istart:iend, jstart:jend, :)
           var_local(:, :, ivar) = sum(var_region, dim=3) / real(ave_days)

           close(iunit)

        end do

        var = reshape(var_local, (/nlon*nlat*var_num/))

     end subroutine get_cfs_field

!----------------------------------------------------------------------
!     Public subroutine to get observed precipitation anomoly percentage
!     in a given date, leading days and average days
!----------------------------------------------------------------------
     subroutine get_prec_anom_direct(input_dir, path_separator, stn_num, year, mon, day, start_year, end_year, fct_mon, fct_day, lead_days, ave_days, prec_anom, ierr)

        implicit none

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
        character(len=*), intent(in)       ::  input_dir
        character(len=*), intent(in)       ::  path_separator
        integer, intent(in)                ::  stn_num
        integer, intent(in)                ::  year, mon, day
        integer, intent(in)                ::  start_year, end_year
        integer, intent(in)                ::  fct_mon, fct_day
        integer, intent(in)                ::  lead_days, ave_days
        real(4), intent(out)               ::  prec_anom(stn_num)
        integer, intent(out)               ::  ierr

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(get_prec_anom_direct)'
        character(256)                     ::  filename
        character(4)                       ::  start_year_str, end_year_str
        character(2)                       ::  fct_mon_str, fct_day_str
        character(2)                       ::  lead_days_str, ave_days_str
        integer                            ::  iunit, status, line
        logical                            ::  lexist, flag
        integer                            ::  year_in, mon_in, day_in
        real(4)                            ::  prec_in(stn_num)

        write( start_year_str, '(i4)' ) start_year
        write( end_year_str, '(i4)' ) end_year
        fct_mon_str = char(fct_mon/10 + 48)//char(mod(fct_mon, 10) + 48)
        fct_day_str = char(fct_day/10 + 48)//char(mod(fct_day, 10) + 48)
        lead_days_str = char(lead_days/10 + 48)//char(mod(lead_days, 10) + 48)
        ave_days_str = char(ave_days/10 + 48)//char(mod(ave_days, 10) + 48)
        ierr = 0
        iunit = 104
        filename = trim(input_dir)//path_separator//start_year_str//'-'//end_year_str//'-'//fct_mon_str//'-'//fct_day_str//'-lead-'//lead_days_str//'-for-'//ave_days_str//'-ano.dat'

        inquire( file = filename, exist = lexist )
        if ( .not. lexist ) then
           write( 6, '(4a)' ) subname, 'ERROR: ', trim(filename) , ' is NOT available, please check whether or not file exist'
           ierr = 1
           return
        end if

        open( unit = iunit, file = filename, form = 'formatted', status = 'old', action = 'read', iostat = status )
        if ( status /= 0 ) then
           write( 6, '(3a)' ) subname, 'ERROR: can NOT rightly open file ', trim(filename)
           ierr = 1
           return
        end if

        line = 1
        flag = .false.
        do while ( .true. )

           read( iunit, *, iostat = status ) year_in, mon_in, day_in, prec_in
           if ( status /= 0 ) then
              write( 6, '(4a, i)' ) subname, 'ERROR: can NOT read precipitation in ', trim(filename), ' for line ', line
              ierr = 1
              return
           end if

           if ( year_in == year .and. mon_in == mon .and. day_in == day ) then

              prec_anom = prec_in
              flag = .true.
              exit

           end if

           line = line + 1

        end do

        close (iunit)

        if ( flag == .false. ) then
           write( 6, '(4a, i4, 1x, i2, 1x, i2)' ) subname, 'ERROR: can NOT find precipitation data in ', trim(filename), ' for ', year, mon, day
           ierr = 1
           return
        end if

     end subroutine get_prec_anom_direct

!----------------------------------------------------------------------
!    Public subroutine to get observed precipitation anmoly percentage
!    in a given date, leading days and average days
!----------------------------------------------------------------------
     subroutine get_prec_anom(input_dir, path_separator, stn_num, year, mon, day, lead_days, ave_days, prec_anom, ierr)

        use date_mod, only :  new_date

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        character(len=*), intent(in)       ::  input_dir
        character(len=*), intent(in)       ::  path_separator
        integer, intent(in)                ::  stn_num
        integer, intent(in)                ::  year, mon, day
        integer, intent(in)                ::  lead_days, ave_days
        real(4), intent(out)               ::  prec_anom(stn_num)
        integer, intent(out)               ::  ierr

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        integer                            ::  start_year, start_mon, start_day
        real(4)                            ::  prec_hist(stn_num, 30)
        real(4)                            ::  prec(stn_num)
        real(4)                            ::  prec_mean(stn_num)
        integer                            ::  iyear

        call new_date(year, mon, day, lead_days, start_year, start_mon, start_day)

        if ( start_mon == 2 .and. start_day == 29 ) then
           start_mon = 3
           start_day = 1
        end if

        do iyear = 1981, 2010
           call get_prec(input_dir, path_separator, stn_num, iyear, start_mon, start_day, ave_days, prec_hist(:, iyear-1980), ierr)
           if ( ierr /= 0 ) then
              return
           end if
        end do

        prec_mean = sum(prec_hist, dim=2) / 30.

        call get_prec(input_dir, path_separator, stn_num, year, start_mon, start_day, ave_days, prec, ierr)

        if ( ierr /= 0 ) then
           return
        end if

        prec_anom = 100. * (prec - prec_mean) / prec_mean

     end subroutine get_prec_anom

!----------------------------------------------------------------------
!    Private subroutine to get observed precipitation in a given date,
!    leading days and average days
!----------------------------------------------------------------------
     subroutine get_prec(input_dir, path_separator, stn_num, year, mon, day, ave_days, prec, ierr)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        character(len=*), intent(in)       ::  input_dir
        character(len=*), intent(in)       ::  path_separator
        integer, intent(in)                ::  stn_num
        integer, intent(in)                ::  year, mon, day
        integer, intent(in)                ::  ave_days
        real(4), intent(out)               ::  prec(stn_num)
        integer, intent(out)               ::  ierr

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(get_prec)'
        character(256)                     ::  filename
        integer                            ::  iunit, status, line
        logical                            ::  lexist, flag
        integer                            ::  year_in, mon_in, day_in
        real(4)                            ::  prec_in(stn_num), prec_days_in(stn_num, ave_days)
        integer                            ::  iday

        ierr = 0
        iunit = 102
        filename = trim(input_dir)//path_separator//'523stn-daily-1980-2013.txt'

        inquire( file = filename, exist = lexist )
        if ( .not. lexist ) then
           write( 6, '(4a)' ) subname, 'ERROR: ', trim(filename) , ' is NOT available, please check whether or not file exist'
           ierr = 1
           return
        end if

        open( unit = iunit, file = filename, form = 'formatted', status = 'old', action = 'read', iostat = status )
        if ( status /= 0 ) then
           write( 6, '(3a)' ) subname, 'ERROR: can NOT rightly open file ', trim(filename)
           ierr = 1
           return
        end if

        line = 1
        flag = .false.
        do while ( .true. )

           read( iunit, *, iostat = status ) year_in, mon_in, day_in, prec_in
           if ( status /= 0 ) then
              write( 6, '(4a, i)' ) subname, 'ERROR: can NOT read precipitation in ', trim(filename), ' for line ', line
              ierr = 1
              return
           end if

           if ( year_in == year .and. mon_in == mon .and. day_in == day ) then

              prec_days_in(:, 1) = prec_in
              flag = .true.
              do iday = 2, ave_days

                 read( iunit, *, iostat = status ) year_in, mon_in, day_in, prec_in
                 if ( status /= 0 ) then
                    write( 6, '(4a, i)' ) subname, 'ERROR: can NOT read precipitation in ', trim(filename), ' for line ', line
                    ierr = 1
                    return
                 end if

                 prec_days_in(:, iday) = prec_in
                 line = line + 1

              end do
              exit

           end if

           line = line + 1

        end do

        close( iunit )

        if ( flag == .false. ) then
           write( 6, '(4a, i4, 1x, i2, 1x, i2)' ) subname, 'ERROR: can NOT find precipitation data in ', trim(filename), ' for ', year, mon, day
           ierr = 1
           return
        end if

        prec = sum(prec_days_in, dim=2) / real(ave_days)

     end subroutine get_prec

  end module inout_mod
