program batch_loglikelihood
  use likelihood_mod, only : readnicerlc, loglik
  use mpi
  implicit none

  ! Parameters
  character(len=*), parameter :: data_dir = 'data/'
  character(len=*), parameter :: file_prefix = 'mcmc_vac'
  character(len=*), parameter :: output_dir = 'output/'
  integer, parameter :: expected_bins = 64
  integer, parameter :: expected_params = 11
  integer, parameter :: nparam_bins = 100
  real(8), parameter :: backgr = 13500.0_8

  ! Variables
  character(len=256), allocatable :: file_list(:)
  integer :: n_files, i, n_to_process, ierr, rank, nprocs, global_best_rank
  character(len=256) :: arg1
  logical :: use_all

  ! Data arrays
  real(8), allocatable :: datalc(:, :), errorlc(:)
  integer :: ndata, ndataact

  ! For per-file results
  real(8), allocatable :: best_model(:), best_params(:)
  real(8) :: best_loglik
  real(8), allocatable :: param_loglik_sum(:,:)

  ! Global best tracking
  real(8) :: global_best_loglik
  real(8), allocatable :: global_best_model(:), global_best_params(:)
  real(8), allocatable :: global_param_loglik_sum(:,:)
  logical :: first

  ! MPI and reduction temporaries
  real(8) :: all_best_loglik
  integer :: best_owner
  real(8), allocatable :: global_param_loglik_sum_all(:,:)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! --- Parse arguments ---
  arg1 = ''
  call get_command_argument(1, arg1)
  use_all = .false.
  n_to_process = 1
  if (trim(arg1) == 'full') then
    use_all = .true.
  else if (len_trim(arg1) > 0) then
    read(arg1, *, iostat=i) n_to_process
    if (i /= 0 .or. n_to_process <= 0) n_to_process = 1
  end if

  ! --- List files (only on rank 0), then broadcast to all ranks ---
  if (rank == 0) then
    call list_files_with_prefix(data_dir, file_prefix, file_list, n_files)
    if (n_files == 0) then
      write(*,*) '[batch_loglikelihood] No files found.'
    end if
  else
    n_files = 0
  end if
  ! Broadcast n_files to all ranks
  call MPI_Bcast(n_files, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (n_files == 0) stop '[batch_loglikelihood] No files found.'
  ! Allocate file_list on all ranks
  if (rank /= 0) allocate(file_list(n_files))
  ! Broadcast file names (as a big character array)
  call bcast_file_list(file_list, n_files, 256, rank, ierr)
  if (use_all) n_to_process = n_files
  if (n_to_process > n_files) n_to_process = n_files
  ! --- Read observed data ---
  ndata = expected_bins
  allocate(datalc(ndata, 2), errorlc(ndata))
  call readnicerlc(datalc, errorlc, ndata, ndataact, backgr)
  if (ndataact <= 0) stop 'No observed data.'

  ! --- Main loop over files ---
  first = .true.
  allocate(global_best_model(ndataact), global_best_params(expected_params))
  allocate(global_param_loglik_sum(expected_params, nparam_bins))
  global_best_loglik = -1.0d99
  global_param_loglik_sum = 0.0d0
  do i = 1, n_to_process
    if (mod(i-1, nprocs) == rank) then
      call process_model_file(trim(data_dir)//trim(file_list(i)), datalc, errorlc, ndata, ndataact, &
           best_loglik, best_model, best_params, param_loglik_sum)
      ! Update local best
      if (first .or. best_loglik > global_best_loglik) then
        global_best_loglik = best_loglik
        global_best_model = best_model
        global_best_params = best_params
        first = .false.
      end if
      ! Accumulate parameter bin loglikelihoods
      global_param_loglik_sum = global_param_loglik_sum + param_loglik_sum
    end if
  end do

  ! --- MPI reduction for best loglikelihood ---
  call MPI_Allreduce(global_best_loglik, all_best_loglik, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  ! Find which rank owns the global best
  if (abs(global_best_loglik - all_best_loglik) < 1d-10) then
    best_owner = rank
  else
    best_owner = -1
  end if
  call MPI_Allreduce(best_owner, global_best_rank, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

  ! Broadcast best params and model from owner to all
  call MPI_Bcast(global_best_params, expected_params, MPI_DOUBLE_PRECISION, global_best_rank, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(global_best_model, ndataact, MPI_DOUBLE_PRECISION, global_best_rank, MPI_COMM_WORLD, ierr)

  ! --- MPI reduction for parameter bin loglikelihoods ---
  if (rank == 0) then
    allocate(global_param_loglik_sum_all(expected_params, nparam_bins))
  end if
  call MPI_Reduce(global_param_loglik_sum, global_param_loglik_sum_all, expected_params*nparam_bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! Only rank 0 prints and writes output
  if (rank == 0) then
    write(*, '(A,ES24.14)') 'Best loglikelihood: ', all_best_loglik
    write(*, '(A)') 'Best parameters:'
    do i = 1, size(global_best_params)
      write(*, '(A,I0,A,ES24.14)') '  p', i, ' = ', global_best_params(i)
    end do
    write(*, '(A)') 'Best model lightcurve:'
    do i = 1, size(global_best_model)
      write(*, '(I8,1X,ES24.14)') i, global_best_model(i)
    end do
    call write_best_lightcurve(global_best_model, 'output/best_model_lightcurve.dat')
    call write_best_params(global_best_params, 'output/best_model_params.dat')
    call write_param_binning(global_param_loglik_sum_all, 'output/param_binning_output.dat')
  end if


  call MPI_Finalize(ierr)

contains

  subroutine bcast_file_list(file_list, n_files, maxlen, rank, ierr)
    character(len=*), allocatable, intent(inout) :: file_list(:)
    integer, intent(in) :: n_files, maxlen, rank
    integer, intent(out) :: ierr
    character(len=maxlen) :: buffer(n_files)
    integer :: i
    if (rank == 0) then
      do i = 1, n_files
        buffer(i) = file_list(i)
      end do
    end if
    call MPI_Bcast(buffer, n_files*maxlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (rank /= 0) then
      do i = 1, n_files
        file_list(i) = buffer(i)
      end do
    end if
  end subroutine bcast_file_list

  subroutine write_best_lightcurve(model, path)
    real(8), intent(in) :: model(:)
    character(len=*), intent(in) :: path
    integer :: unit, ios, i
    open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    do i = 1, size(model)
      write(unit, '(F10.6,1X,ES24.14)') real(i-1,8)/real(expected_bins,8), model(i)
    end do
    close(unit)
  end subroutine write_best_lightcurve

  subroutine write_best_params(params, path)
    real(8), intent(in) :: params(:)
    character(len=*), intent(in) :: path
    integer :: unit, ios, i
    open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    do i = 1, size(params)
      write(unit, '(A,I0,A,ES24.14)') 'p', i, ' ', params(i)
    end do
    close(unit)
  end subroutine write_best_params

  subroutine write_param_binning(param_loglik_sum, path)
    real(8), intent(in) :: param_loglik_sum(:,:)
    character(len=*), intent(in) :: path
    integer :: unit, ios, p, bin
    open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(unit, '(A)') '# Parameter binning output: param, bin_idx, loglik_sum'
    do p = 1, size(param_loglik_sum,1)
      do bin = 1, size(param_loglik_sum,2)
        write(unit, '(I3,1X,I3,1X,ES24.14)') p, bin, param_loglik_sum(p,bin)
      end do
    end do
    close(unit)
  end subroutine write_param_binning

  subroutine list_files_with_prefix(dir, prefix, files, n)
    character(len=*), intent(in) :: dir, prefix
    character(len=256), allocatable, intent(out) :: files(:)
    integer, intent(out) :: n
    character(len=4096) :: cmd, line
    integer :: ios, unit, count
    cmd = 'ls -1 ' // trim(dir) // trim(prefix) // '* > filelist.tmp'
    call execute_command_line(cmd)
    open(newunit=unit, file='filelist.tmp', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      n = 0
      return
    end if
    count = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      count = count + 1
    end do
    rewind(unit)
    allocate(files(count))
    n = count
    count = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      count = count + 1
      ! Remove directory prefix if present
      if (index(trim(line), trim(dir)) == 1) then
        files(count) = adjustl(trim(line(len_trim(dir)+1:)))
      else
        files(count) = adjustl(trim(line))
      end if
    end do
    close(unit)
    call execute_command_line('rm -f filelist.tmp')
  end subroutine list_files_with_prefix

  subroutine process_model_file(model_path, datalc, errorlc, ndata, ndataact, &
       best_loglik, best_model, best_params, param_loglik_sum)
    character(len=*), intent(in) :: model_path
    integer, intent(in) :: ndata, ndataact
    real(8), intent(in) :: datalc(ndata, 2), errorlc(ndata)
    real(8), allocatable, intent(out) :: best_model(:), best_params(:)
    real(8), intent(out) :: best_loglik
    real(8), allocatable, intent(out) :: param_loglik_sum(:,:)
    integer, parameter :: expected_total_cols = 76
    integer :: unit, ios, nparams, nmodels, case_idx, bin_idx, p
    real(8), allocatable :: model(:), params(:)
    real(8) :: rowvals(expected_total_cols), current_ll
    real(8), allocatable :: param_minvals(:), param_maxvals(:)
    integer :: nparam_bins
    nparam_bins = 100
    nparams = 11
    allocate(model(ndataact), best_model(ndataact), params(nparams), best_params(nparams))
    allocate(param_loglik_sum(nparams, nparam_bins))
    allocate(param_minvals(nparams), param_maxvals(nparams))
    param_loglik_sum = 0.0d0
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
    best_loglik = -1.0d99
    nmodels = 0
    open(newunit=unit, file=model_path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) '[process_model_file] Could not open model file:', trim(model_path), ' iostat=', ios
      return
    end if
    do
      rowvals = 0.0d0
      read(unit, *, iostat=ios) rowvals
      if (ios /= 0) exit
      params = rowvals(1:nparams)
      model = rowvals(13:76)
      current_ll = loglik(datalc, errorlc, model, ndata, ndataact)
      nmodels = nmodels + 1
      if (current_ll > best_loglik) then
        best_loglik = current_ll
        best_model = model
        best_params = params
      end if
      do p = 1, nparams
        if (param_maxvals(p) > param_minvals(p)) then
          bin_idx = int((params(p) - param_minvals(p)) / (param_maxvals(p) - param_minvals(p)) * (nparam_bins - 1)) + 1
          if (bin_idx < 1) bin_idx = 1
          if (bin_idx > nparam_bins) bin_idx = nparam_bins
          param_loglik_sum(p, bin_idx) = param_loglik_sum(p, bin_idx) + exp(current_ll + 2100.0d0)
        end if
      end do
    end do
    close(unit)
  end subroutine process_model_file

  subroutine write_batch_outputs(filename, best_loglik, best_params, best_model, param_loglik_sum)
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: best_loglik
    real(8), intent(in) :: best_params(:), best_model(:)
    real(8), intent(in) :: param_loglik_sum(:,:)
    integer :: unit, ios, p, bin, nparam_bins, nparams
    character(len=256) :: outpath
    nparams = size(best_params)
    nparam_bins = size(param_loglik_sum, 2)
    ! Remove any directory prefix from filename for output
    if (index(trim(filename), '/') > 0) then
      outpath = trim(output_dir)//trim(filename(index(trim(filename),'/')+1:))//'_summary.dat'
    else
      outpath = trim(output_dir)//trim(filename)//'_summary.dat'
    end if
    open(newunit=unit, file=outpath, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) '[write_batch_outputs] Could not open output file:', trim(outpath), ' iostat=', ios
      return
    end if
    write(*,*) '[write_batch_outputs] Writing output file:', trim(outpath)
    write(unit, '(A,F18.6)') 'Best loglikelihood: ', best_loglik
    write(unit, '(A)') 'Best parameters:'
    do p = 1, nparams
      write(unit, '(A,I0,A,ES24.14)') '  p', p, ' = ', best_params(p)
    end do
    write(unit, '(A)') 'Best model lightcurve:'
    do p = 1, size(best_model)
      write(unit, '(I8,1X,ES24.14)') p, best_model(p)
    end do
    write(unit, '(A)') 'Parameter binning (param, bin, value):'
    do p = 1, nparams
      do bin = 1, nparam_bins
        write(unit, '(I3,1X,I3,1X,E32.16)') p, bin, param_loglik_sum(p, bin)
      end do
    end do
    close(unit)
  end subroutine write_batch_outputs

end program batch_loglikelihood

