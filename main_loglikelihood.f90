program loglikelihood_main
  use likelihood_mod, only : readnicerlc, evaluate_model_file, preview_model_file_head, write_lightcurve_files
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none

  character(len=*), parameter :: data_path = 'nicer_profile.dat'
  character(len=*), parameter :: default_model_path = '/home7/jraynau1/nobackup/loglikelihood/640m_parameters_and_phase_amplitudes.dat'
  real(8), parameter :: backgr = 13500.0_8
  integer, parameter :: expected_bins = 64
  integer, parameter :: default_cases_to_test = 1
  character(len=*), parameter :: nicer_out_path = 'output/nicer_lightcurve.dat'
  character(len=*), parameter :: best_model_out_path = 'output/best_model_lightcurve.dat'
  character(len=*), parameter :: first_model_out_path = 'output/first_model_lightcurve.dat'

  integer :: ndata, ndataact, nparams, nmodels, i, cases_to_test, parse_ios
  real(8), allocatable :: datalc(:, :), errorlc(:), best_model(:), first_model(:), best_params(:)
  real(8) :: best_ll
  real(8) :: t0, t1
  character(len=1024) :: arg1, arg2, arg3, model_path
  logical :: full_run

  arg1 = ''
  arg2 = ''
  arg3 = ''
  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)

  full_run = .true.
  model_path = default_model_path
  cases_to_test = default_cases_to_test

  call log_msg('[main] Program started')

  if (trim(arg1) == '--preview') then
    full_run = .false.
    if (len_trim(arg2) > 0) model_path = arg2
  else if (trim(arg1) == '--full') then
    full_run = .true.
    if (len_trim(arg2) > 0) model_path = arg2
    if (len_trim(arg3) > 0) then
      read(arg3, *, iostat=parse_ios) cases_to_test
      if (parse_ios /= 0 .or. cases_to_test <= 0) cases_to_test = default_cases_to_test
    end if
  else
    if (len_trim(arg1) > 0) model_path = arg1
    if (len_trim(arg2) > 0) then
      read(arg2, *, iostat=parse_ios) cases_to_test
      if (parse_ios /= 0 .or. cases_to_test <= 0) cases_to_test = default_cases_to_test
    end if
  end if

  if (.not. full_run) then
    call log_msg('[main] Preview mode: printing only the first 75 values')
    call preview_model_file_head(trim(model_path))
    call log_msg("[main] Full mode evaluates requested fixed-size cases")
    stop
  end if

  ndata = expected_bins

  allocate(datalc(ndata, 2), errorlc(ndata))
  call log_msg('[main] Reading observed data + errors')
  call cpu_time(t0)
  call readnicerlc(datalc, errorlc, ndata, ndataact, backgr)
  call cpu_time(t1)
  write(*, '(A,F10.3,A)') '[main] readnicerlc done in ', t1 - t0, ' s'

  if (ndataact <= 0) then
    write(*, '(A)') 'No rows were read from '//trim(data_path)
    stop 1
  end if

  call log_msg('[main] Evaluating model file')
  write(output_unit, '(A,I0)') '[main] Cases requested: ', cases_to_test
  call flush(output_unit)
  call cpu_time(t0)
  call evaluate_model_file(trim(model_path), datalc, errorlc, ndata, ndataact, &
       best_ll, best_model, first_model, best_params, nparams, nmodels, cases_to_test)
  call cpu_time(t1)
  write(*, '(A,F10.3,A)') '[main] evaluate_model_file done in ', t1 - t0, ' s'

  call write_lightcurve_files(datalc, best_model, first_model, ndata, ndataact, &
       nicer_out_path, best_model_out_path, first_model_out_path)
  write(output_unit, '(A)') '[main] NICER lightcurve file path: '//trim(nicer_out_path)
  write(output_unit, '(A)') '[main] Best model lightcurve file path: '//trim(best_model_out_path)
  write(output_unit, '(A)') '[main] First model lightcurve file path: '//trim(first_model_out_path)
  call flush(output_unit)

  call log_msg('Best model search complete.')
  write(*, '(A)') 'Model file: '//trim(model_path)
  write(*, '(A,I0)') 'Observed rows: ', ndataact
  write(*, '(A,I0)') 'Case limit: ', cases_to_test
  write(*, '(A,I0)') 'Models evaluated: ', nmodels
  write(*, '(A,F18.6)') 'Best log-likelihood: ', best_ll

  if (nparams > 0) then
    write(*, '(A)') 'Best parameter set:'
    do i = 1, nparams
      write(*, '(A,I0,A,ES24.14)') '  p', i, ' = ', best_params(i)
    end do
  end if

  call log_msg('Best lightcurve values:')
  do i = 1, ndataact
    write(*, '(I8,1X,ES24.14)') i, best_model(i)
  end do

contains

  subroutine log_msg(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(output_unit, '(A)') trim(msg)
    call flush(output_unit)
  end subroutine log_msg

end program loglikelihood_main

