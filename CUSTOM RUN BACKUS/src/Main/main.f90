program main
    use problem_module, only : problem_t
    use datafile_module, only : datafile_t
    use mpi
    use omp_lib
    implicit none

    type(problem_t) , allocatable    :: prob
    type(datafile_t) :: df_obj

    integer :: num_args, ix
    integer :: ierr, rank

    character(:), allocatable :: arg
    character(:), allocatable :: env_datafile
    integer :: arglen, stat
    integer :: env_len

    call get_environment_variable("SCALSALE_DATAFILE", length=env_len, status=stat)
    if (env_len > 0 .and. stat == 0) then
        allocate(character(env_len) :: env_datafile)
        call get_environment_variable("SCALSALE_DATAFILE", value=env_datafile, status=stat)
    end if

    call get_command_argument(number=1, length=arglen)  ! Assume for simplicity success
    if (arglen == 0) then
        if (allocated(env_datafile)) then
            df_obj = datafile_t(env_datafile)
        else
            df_obj = datafile_t("../Datafiles/datafile.json")
        end if
        else
            allocate (character(arglen) :: arg)
            call get_command_argument(number=1, value=arg, status=stat)
            df_obj = datafile_t(arg)
    end if

    call MPI_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    allocate(prob)



    prob = problem_t(df_obj)
    call prob%Start_calculation()
    call MPI_FINALIZE(ierr)

end program main
