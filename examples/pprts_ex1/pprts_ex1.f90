module m_pprts_ex1
    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi
    use m_pprts, only : init_pprts, set_optical_properties, t_solver, t_solver_3_6, t_solver_8_10, &
      solve_pprts, pprts_get_result, set_angles
    use m_tenstream_options, only: read_commandline_options

    use mpi, only : MPI_COMM_WORLD

  implicit none

  contains
subroutine pprts_ex1()
    integer(iintegers),parameter :: nxp=9,nyp=9,nv=2
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=270, theta0=20
    real(ireals),parameter :: albedo=0., dz=100
    real(ireals),parameter :: incSolar = 1364
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver

    character(len=80) :: arg

    dz1d = dz

    call get_command_argument(1, arg)

    select case(arg)
      case ('-solver_8_10')
        allocate(t_solver_8_10::solver)
      case('-solver_3_6')
        allocate(t_solver_3_6::solver)
      case default
        print *,'error, have to provide solver type as argument, e.g.'
        print *,'-solver_8_10'
        print *,'-solver_3_6'
        stop
    end select

    call init_pprts(mpi_comm_world, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = .1_ireals/(dz*nv)
    ksca = 1e-3_ireals/(dz*nv)
    g    = zero

    kabs(nv/2,nxp/2,1:nyp) = 1/dz
    ksca(nv/2,nxp/2,1:nyp) = 1/dz
    g   (nv/2,nxp/2,1:nyp) = .9

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call set_angles(solver, 0._ireals, theta0)

    call solve_pprts(solver, incSolar)


    allocate(fdir (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdn  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fup  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdiv (solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    call pprts_get_result(solver, fdir,fdn,fup,fdiv)

    print *,'edir', fdir(:, nxp/2,nyp/2)
    print *,'edn:', fdn(:, nxp/2,nyp/2)
    print *,'eup:', fup(:, nxp/2,nyp/2)
end subroutine

!subroutine example(solver)
!
!end subroutine
end module


program main
  use m_pprts_ex1
  integer(mpiint) :: myid, ierr

  call pprts_ex1()
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  if(myid.eq.0) then
    print *,''
    print *,''
    print *,'Call this example e.g. with options: -show_edir hdf5:edir.h5'
    print *,'and plot results with python:'
    print *,'import h5py as H; h=H.File("edir.h5","r"); edir = h["edir0"][:]'
    print *,'imshow(edir[0,:,:,0].T,interpolation="nearest");' ! has dimension nyp,nxp,nzp,8streams
    print *,'colorbar(); savefig("edir_x0.pdf")'
  endif

end program