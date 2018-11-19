!***********************************************************************
!  DESCRIPTION:
!       Module to load some necessary parameters to control the whole
!       behavior of dynamical similarity method
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!
!***********************************************************************
  module params_mod

     implicit none

!----------------------------------------------------------------------
!    Parameter definitions for namelist file
!
!    fct_date        = the date when CFS v2 begin to forecast
!    start_year      = start year to search for similar historical forecast
!    end_year        = end year to search for similar historical forecast
!    lead_days       = leading days for dynamical similarity prediction
!    ave_days        = period length for avarage
!    domain_west     = west edge of interested area
!    domain_east     = east edge of interested area
!    domain_south    = south edge of interested area
!    domain_north    = north edge of interested area
!    fct_num         = selected forecast number
!    var_num         = varible number
!    var_names       = varible names
!    input_dir_curt  = directory of current CFSv2 forecast
!    input_dir_hist  = directory of historical CFSv2 forecast
!    input_dir_prec  = directory of historical CFSv2 forecast
!    output_dir      = directory for output
!
!----------------------------------------------------------------------
     integer, parameter                 ::  max_num = 10
     integer                            ::  resample_num
     real(4)                            ::  alpha
     integer                            ::  fct_mon, fct_day
     integer                            ::  lead_days, ave_days
     integer                            ::  domain_west, domain_east, domain_south, domain_north
     integer                            ::  start_year, end_year
     integer                            ::  samp_interval
     integer                            ::  samp_sgl_len
     integer                            ::  anal_fct_num ! selected forecast number
     integer                            ::  var_num ! for further research purpose
     character(256)                     ::  var_names(max_num)
     character(256)                     ::  input_dir_curt, input_dir_hist, input_dir_prec, output_dir
     character(1)                       ::  path_separator
     integer                            ::  debug_level

     public          ::  init_params   !initialize the parameters listed above

     contains

!----------------------------------------------------------------------
!    Public subroutine to initialize control parameters
!----------------------------------------------------------------------
     subroutine init_params(filename, ierr)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        character(len=*), intent(in)       ::  filename   ! namelist file name
        integer, intent(out)               ::  ierr       ! error identifier

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(init_params)'
        integer                            ::  iunit = 255     ! file unit

        namelist /params/ resample_num, alpha, fct_mon, fct_day, lead_days, ave_days, domain_west, domain_east, domain_south, domain_north, &
                          start_year, end_year, samp_interval, samp_sgl_len, anal_fct_num, var_num, var_names, input_dir_curt, input_dir_hist, &
                          input_dir_prec, output_dir, path_separator, debug_level

!------------------------------------------------------------------
!    Open the namelist file and obtain parameter's value
!------------------------------------------------------------------
        open( unit = iunit, file = filename, status = 'old', action = 'read', err = 100 )

        read( iunit, params, err = 200 )

        ierr = 0
        return

100     write( 6, '(3a)' ) subname, 'can NOT read namelist file ', trim(filename)
        ierr = 1
        return

200     write( 6, '(3a)' ) subname, 'can NOT read namelist stanza: params  ',  trim(filename)
        ierr = 1
        close( iunit )

     end subroutine init_params

  end module params_mod
