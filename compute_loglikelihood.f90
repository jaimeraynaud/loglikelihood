module likelihood_mod
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none

contains

  subroutine readnicerlc(datalc, errorlc, ndata, ndataact, backgr)
    implicit none
    integer, intent(in) :: ndata
    integer, intent(out) :: ndataact
    real(8), intent(out) :: datalc(ndata, 2), errorlc(ndata)
    real(8), intent(in) :: backgr
    real(8) :: vmx, x, y
    integer :: m, i, ios
    real(8) :: t0, t1

    call cpu_time(t0)
    write(*, '(A)') '[readnicerlc] Opening nicer_profile.dat'

    open(21, file='nicer_profile.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      ndataact = 0
      write(*, '(A,I0)') '[readnicerlc] Open failed. iostat=', ios
      return
    end if

    write(*, '(A)') '[readnicerlc] Reading observed lightcurve rows'
    m = 1
    do
      read(21, *, iostat=ios) x, y
      if (ios /= 0) exit
      if (m > ndata) exit
      datalc(m, 1) = x
      datalc(m, 2) = y - backgr
      m = m + 1
    end do
    close(21)

    ndataact = m - 1
    write(*, '(A,I0)') '[readnicerlc] Rows read: ', ndataact
    if (ndataact <= 0) return

    vmx = 0.0d0
    do i = 1, ndataact
      if (datalc(i, 2) > vmx) vmx = datalc(i, 2)
    end do

    do i = 1, ndataact
      datalc(i, 2) = datalc(i, 2)
      errorlc(i) = sqrt(datalc(i, 2) + 2.0d0 * backgr)
    end do

    call cpu_time(t1)
    write(*, '(A,F10.3,A)') '[readnicerlc] Done in ', t1 - t0, ' s'
  end subroutine readnicerlc

  real(8) function loglik(datalc, errorlc, modellc, ndata, ndataact)
    implicit none
    integer, intent(in) :: ndata, ndataact
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real(8), intent(in) :: datalc(ndata, 2), errorlc(ndata), modellc(ndataact)
    real(8) :: err2
    integer :: i

    loglik = 0.0d0
    do i = 1, ndataact
      err2 = errorlc(i) ** 2
      loglik = loglik + (datalc(i, 2) - modellc(i)) ** 2 / err2 + log(2.0d0 * pi * err2)
    end do
    loglik = -0.5d0 * loglik
  end function loglik

  subroutine evaluate_model_file(model_path, datalc, errorlc, ndata, ndataact, best_loglik, &
       best_model, first_model, best_params, nparams, nmodels, max_models)
    implicit none
    character(len=*), intent(in) :: model_path
    integer, intent(in) :: ndata, ndataact
    real(8), intent(in) :: datalc(ndata, 2), errorlc(ndata)
    real(8), intent(out) :: best_loglik
    real(8), allocatable, intent(out) :: best_model(:), first_model(:), best_params(:)
    integer, intent(out) :: nparams, nmodels
    integer, intent(in), optional :: max_models

    integer, parameter :: expected_total_cols = 76
    integer, parameter :: expected_params = 11
    integer, parameter :: expected_bins = 64
    integer :: unit, ios, target_cases, case_idx, ii
    real(8), allocatable :: model(:), params(:)
    real(8) :: rowvals(expected_total_cols)
    real(8) :: current_ll, file_ll, rel_error, rel_error_sum, epsilon
    real(8) :: best_file_ll
    real(8) :: t_open0, t_open1
    ! --- Parameter binning setup ---
    integer, parameter :: nparam_bins = 100
    real(8), allocatable :: param_bin_min(:), param_bin_max(:)
    real(8), allocatable :: param_minvals(:), param_maxvals(:)
    real(8), allocatable :: param_loglik_sum(:,:)
    integer :: bin_idx, p, total_cases, bin, fout
    character(len=256) :: param_binning_file

    call cpu_time(t_open0)
    write(*, '(A)') '[evaluate_model_file] Opening model file'
    open(newunit=unit, file=model_path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') 'Could not open model file: '//trim(model_path)
      stop 1
    end if
    call cpu_time(t_open1)
    write(*, '(A,F10.3,A)') '[evaluate_model_ile] Open done in ', t_open1 - t_open0, ' s'
    write(*, '(A)') '[evaluate_model_file] Reading cases as fixed blocks of 76 values (11 params + 64 bins)'
    call flush(output_unit)

    if (ndataact /= expected_bins) then
      close(unit)
      write(*, '(A,I0,A,I0,A)') '[evaluate_model_file] Expected ', expected_bins, ' observed bins, found ', ndataact, '.'
      stop 1
    end if

    ! Process all cases in the file, not just a fixed number
    nparams = expected_params
    allocate(model(ndataact), best_model(ndataact), first_model(ndataact), params(nparams), best_params(nparams))
    nmodels = 0
    rel_error_sum = 0.0d0
    epsilon = 1.0d-12
    best_file_ll = 0.0d0
    best_loglik = -1.0d99
    ! Now best_loglik is properly initialized for max search

    ! Print the first 76 values for inspection
    rowvals = 0.0d0
    read(unit, *, iostat=ios) rowvals
    if (ios /= 0) then
      write(*, '(A)') '[evauate_model_file] Could not read the first case for printing.'
      close(unit)
      stop 1
    end if
    write(*, '(A)') '[evaluate_model_file] First case (76 values):'
    write(*, '(1X,ES16.8)') (rowvals(ii), ii=1, expected_total_cols)

    ! Now process the first case as usual
    params = rowvals(1:11)
    file_ll = rowvals(12)
    model = rowvals(13:76)
    current_ll = loglik(datalc, errorlc, model, ndata, ndataact)
    rel_error = abs(current_ll - file_ll) / max(abs(file_ll), epsilon)
    rel_error_sum = rel_error_sum + rel_error
    first_model = model
    nmodels = 1
    if (current_ll > best_loglik) then
      best_loglik = current_ll
      best_model = model
      best_params = params
      best_file_ll = file_ll
    end if

    ! Print the computed and file loglikelihood for the first 11 cases for inspection
    do case_idx = 1, 11
      rowvals = 0.0d0
      read(unit, *, iostat=ios) rowvals
      if (ios /= 0) exit
      params = rowvals(1:11)
      file_ll = rowvals(12)
      model = rowvals(13:76)
      current_ll = loglik(datalc, errorlc, model, ndata, ndataact)
      write(*, '(A,I0,A,F18.6,A,F18.6)') '[evaluate_model_file] Case ', case_idx, ': computed loglik = ', current_ll, ', file loglik = ', file_ll
    end do

    ! Continue with the rest of the file
    do
      rowvals = 0.0d0
      read(unit, *, iostat=ios) rowvals
      if (ios /= 0) exit
      params = rowvals(1:11)
      file_ll = rowvals(12)
      model = rowvals(13:76)
      current_ll = loglik(datalc, errorlc, model, ndata, ndataact)
      rel_error = abs(current_ll - file_ll) / max(abs(file_ll), epsilon)
      rel_error_sum = rel_error_sum + rel_error
      nmodels = nmodels + 1

      if (current_ll > best_loglik) then
        best_loglik = current_ll
        best_model = model
        best_params = params
        best_file_ll = file_ll
      end if
    end do

    ! close(unit) -- moved to after all file operations

    if (nmodels <= 0) then
      write(*, '(A)') '[evaluate_model_file] Could not read any complete 76-value case from file start'
      stop 1
    end if

    if (present(max_models)) then
      write(*, '(A,I0)') '[evaluate_model_file] Requested cases: ', max_models
    end if
    write(*, '(A,I0)') '[evaluate_model_file] Cases evaluated: ', nmodels
    write(*, '(A,F18.6)') '[evaluate_model_file] Best loglik (computed) = ', best_loglik
    write(*, '(A,F18.6)') '[evaluate_model_file] Best loglik (from file) = ', best_file_ll
    if (nmodels > 0) then
      write(*, '(A,ES12.5)') '[evaluate_model_file] Relative error for best case: ', abs(best_loglik - best_file_ll) / max(abs(best_file_ll), epsilon)
    end if
    call flush(output_unit)

    allocate(param_loglik_sum(expected_params, nparam_bins))
                    allocate(param_minvals(expected_params), param_maxvals(expected_params))
                    ! Set hardcoded min/max values for each parameter
                    param_minvals = 0.0d0
                    param_maxvals = 0.0d0
                    param_maxvals(1) = 0.7d0
                    param_maxvals(2) = 0.7d0
                    param_maxvals(3) = 0.7d0
                    param_maxvals(4) = 3.141592653589793d0
                    param_maxvals(5) = 6.283185307179586d0
                    param_maxvals(6) = 0.7d0
                    param_maxvals(7) = 0.7d0
                    param_maxvals(8) = 0.7d0
                    param_maxvals(9) = 3.141592653589793d0
                    param_maxvals(10) = 6.283185307179586d0
                    param_maxvals(11) = 12.0d0
    param_loglik_sum = 0.0d0

    ! --- Bin loglikelihoods using hardcoded min/max values ---
    rewind(unit)
    do
      rowvals = 0.0d0
      read(unit, *, iostat=ios) rowvals
      if (ios /= 0) exit
      params = rowvals(1:expected_params)
      model = rowvals(expected_params+2:expected_total_cols)
      current_ll = loglik(datalc, errorlc, model, ndata, ndataact)
      do p = 1, expected_params
        if (param_maxvals(p) > param_minvals(p)) then
          bin_idx = int((params(p) - param_minvals(p)) / (param_maxvals(p) - param_minvals(p)) * (nparam_bins - 1)) + 1
          if (bin_idx < 1) bin_idx = 1
          if (bin_idx > nparam_bins) bin_idx = nparam_bins
          !Old version:
          !param_loglik_sum(p, bin_idx) = param_loglik_sum(p, bin_idx) + current_ll
          !New version:
          param_loglik_sum(p, bin_idx) = param_loglik_sum(p, bin_idx) + exp(current_ll + 2100.0d0)
        end if
      end do
    end do

    ! --- Output binned loglikelihoods ---
    do p = 1, expected_params
      write(*, '(A,I0,A,2ES16.8)') '[param_binning] Parameter ', p, ' min/max: ', param_minvals(p), param_maxvals(p)
      write(*, '(A,I0)') '[param_binning] Bin loglikelihood sums for parameter ', p
      write(*, '(100ES16.8)') (param_loglik_sum(p, bin_idx), bin_idx=1,nparam_bins)
    end do

    ! --- Write parameter binning results to file for Python plotting ---
    close(unit)
    param_binning_file = 'output/param_binning_output.dat'
    open(newunit=fout, file=param_binning_file, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') '[param_binning] Could not open output file.'
    else
      write(fout, '(A)') '# Parameter bining output: param, bin_idx, bin_left, bin_right, loglik_sum'
      do p = 1, expected_params
        do bin = 1, nparam_bins
          write(fout, '(I3,1X,I3,1X,E32.16,1X,E32.16,1X,E32.16,1X,E32.16)') p, bin, &
            param_minvals(p) + (bin-1)*(param_maxvals(p)-param_minvals(p))/nparam_bins, &
            param_minvals(p) + bin*(param_maxvals(p)-param_minvals(p))/nparam_bins, &
            param_loglik_sum(p, bin), &
            0.0d0
        end do
      end do
      close(fout)
      write(*, '(A)') '[param_binning] Wrote parameter binning results to param_binning_output.dat'
    end if

  end subroutine evaluate_model_file

  subroutine write_lightcurve_files(datalc, best_model_lc, first_model_lc, ndata, ndataact, &
       nicer_out_path, best_model_out_path, first_model_out_path)
    implicit none
    integer, intent(in) :: ndata, ndataact
    real(8), intent(in) :: datalc(ndata, 2), best_model_lc(ndataact), first_model_lc(ndataact)
    character(len=*), intent(in) :: nicer_out_path, best_model_out_path, first_model_out_path

    integer :: unit_obs, unit_best, unit_first, ios, i
    logical :: exists_obs, exists_best, exists_first
    character(len=2048) :: cwd, full_obs_path, full_best_path, full_first_path

    open(newunit=unit_obs, file=nicer_out_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') '[write_lightcurve_files] Could not open output: '//trim(nicer_out_path)
      return
    end if

    open(newunit=unit_best, file=best_model_out_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      close(unit_obs)
      write(*, '(A)') '[write_lightcurve_files] Could not open output: '//trim(best_model_out_path)
      return
    end if

    open(newunit=unit_first, file=first_model_out_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      close(unit_obs)
      close(unit_best)
      write(*, '(A)') '[write_lightcurve_files] Could not open output: '//trim(first_model_out_path)
      return
    end if

    do i = 1, ndataact
      write(unit_obs, '(ES24.14,1X,ES24.14)') datalc(i, 1), datalc(i, 2)
      write(unit_best, '(ES24.14,1X,ES24.14)') datalc(i, 1), best_model_lc(i)
      write(unit_first, '(ES24.14,1X,ES24.14)') datalc(i, 1), first_model_lc(i)
    end do

    close(unit_obs)
    close(unit_best)
    close(unit_first)

    call get_environment_variable('PWD', cwd)
    if (len_trim(cwd) > 0) then
      full_obs_path = trim(cwd)//'/'//trim(nicer_out_path)
      full_best_path = trim(cwd)//'/'//trim(best_model_out_path)
      full_first_path = trim(cwd)//'/'//trim(first_model_out_path)
    else
      full_obs_path = trim(nicer_out_path)
      full_best_path = trim(best_model_out_path)
      full_first_path = trim(first_model_out_path)
    end if

    inquire(file=trim(nicer_out_path), exist=exists_obs)
    inquire(file=trim(best_model_out_path), exist=exists_best)
    inquire(file=trim(first_model_out_path), exist=exists_first)

    if (exists_obs) then
      write(*, '(A)') '[write_lightcurve_files] Wrote observed lightcurve: '//trim(full_obs_path)
    else
      write(*, '(A)') '[write_lightcurve_files] ERROR: observed lightcurve file not found after write: '//trim(full_obs_path)
    end if

    if (exists_best) then
      write(*, '(A)') '[write_lightcurve_files] Wrote best lightcurve: '//trim(full_best_path)
    else
      write(*, '(A)') '[write_lightcurve_files] ERROR: best lightcurve file not found after write: '//trim(full_best_path)
    end if

    if (exists_first) then
      write(*, '(A)') '[write_lightcurve_files] Wrote first lightcurve: '//trim(full_first_path)
    else
      write(*, '(A)') '[write_lightcurve_files] ERROR: first lightcurve file not found after write: '//trim(full_first_path)
    end if
  end subroutine write_lightcurve_files

  subroutine preview_model_file_head(model_path)
    implicit none
    character(len=*), intent(in) :: model_path

    integer, parameter :: values_per_row = 75
    integer, parameter :: total_needed = 75
    integer :: unit, ios, parse_ios, ntok, take, collected
    integer(8) :: lines_scanned
    character(len=1048576) :: line
    real(8) :: first_values(total_needed)
    real(8), allocatable :: line_values(:)

    write(*, '(A)') '[preview_model_file_head] Opening model file'
    open(newunit=unit, file=model_path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') '[preview_model_file_head] Could not open model file: '//trim(model_path)
      return
    end if

    first_values = 0.0d0
    collected = 0
    lines_scanned = 0_8
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      lines_scanned = lines_scanned + 1_8
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ntok = count_tokens(line)
      if (ntok <= 0) cycle

      allocate(line_values(ntok))
      line_values = 0.0d0
      read(line, *, iostat=parse_ios) line_values
      if (parse_ios /= 0) then
        deallocate(line_values)
        cycle
      end if

      take = min(ntok, total_needed - collected)
      first_values(collected + 1:collected + take) = line_values(1:take)
      collected = collected + take
      deallocate(line_values)

      if (collected >= total_needed) exit
    end do

    if (collected <= 0) then
      write(*, '(A,I0)') '[preview_model_file_head] No usable numeric values found. Lines scanned: ', lines_scanned
      close(unit)
      return
    end if

    write(*, '(A,I0)') '[preview_model_file_head] Lines scanned: ', lines_scanned
    write(*, '(A,I0,A,I0)') '[preview_model_file_head] Values collected: ', collected, ' / ', total_needed

    if (collected >= values_per_row) then
      write(*, '(A)') '[preview_model_file_head] First 75 values:'
      write(*, '(*(ES16.8,1X))') (first_values(take), take = 1, values_per_row)
    else
      write(*, '(A)') '[preview_model_file_head] First-row values (partial):'
      write(*, '(*(ES16.8,1X))') (first_values(take), take = 1, collected)
    end if

    close(unit)
  end subroutine preview_model_file_head

  subroutine print_rows_for_fixed_columns(model_path, fixed_cols)
    implicit none
    character(len=*), intent(in) :: model_path
    integer, intent(in) :: fixed_cols

    integer :: unit, ios, ntok
    integer(8) :: total_tokens, lines_scanned, full_rows, remainder
    character(len=1048576) :: line

    if (fixed_cols <= 0) then
      write(*, '(A)') '[print_rows_for_fixed_columns] fixed_cols must be > 0'
      return
    end if

    write(*, '(A)') '[print_rows_for_fixed_columns] Opening model file'
    open(newunit=unit, file=model_path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') '[print_rows_for_fixed_columns] Could not open model file: '//trim(model_path)
      return
    end if

    total_tokens = 0_8
    lines_scanned = 0_8
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      lines_scanned = lines_scanned + 1_8
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle
      ntok = count_tokens(line)
      total_tokens = total_tokens + int(ntok, kind=8)
    end do
    close(unit)

    full_rows = total_tokens / int(fixed_cols, kind=8)
    remainder = mod(total_tokens, int(fixed_cols, kind=8))

    write(*, '(A,I0)') '[print_rows_for_fixed_columns] Lines scanned: ', lines_scanned
    write(*, '(A,I0)') '[print_rows_for_fixed_columns] Total values: ', total_tokens
    write(*, '(A,I0,A,I0)') '[print_rows_for_fixed_columns] Rows at ', fixed_cols, ' values/row: ', full_rows
    if (remainder /= 0_8) then
      write(*, '(A,I0)') '[print_rows_for_fixed_columns] Leftover values (incomplete row): ', remainder
    end if
  end subroutine print_rows_for_fixed_columns

  integer function count_tokens(line)
    implicit none
    character(len=*), intent(in) :: line
    integer :: i
    logical :: in_token

    count_tokens = 0
    in_token = .false.

    do i = 1, len_trim(line)
      if (line(i:i) == ' ' .or. line(i:i) == char(9) .or. line(i:i) == ',') then
        in_token = .false.
      else
        if (.not. in_token) count_tokens = count_tokens + 1
        in_token = .true.
      end if
    end do
  end function count_tokens

end module likelihood_mod
