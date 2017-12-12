module test_pprts_symmetry

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, zero, one, pi, mpiint

#include "petsc/finclude/petsc.h"
  use petsc

  use m_pprts, only : init_pprts, set_optical_properties, t_solver_3_6, t_solver_8_10, &
    solve_pprts, set_angles, pprts_get_result, destroy_pprts
  use m_tenstream_options, only: read_commandline_options

  use m_optprop, only: t_optprop, t_optprop_8_10, t_optprop_3_6
  use pfunit_mod

  implicit none

contains

  !@test(npes = [2])
  subroutine test_pprts_symmetry_ex2(this)

    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: numnodes, comm
    integer(iintegers) :: myid
    integer(mpiint) :: ierr

    real(ireals),allocatable  :: dir2dir(:), dir2diff(:)
    integer(iintegers) :: i

    class(t_optprop),allocatable :: OPP

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    allocate (t_optprop_8_10 :: OPP)

    call this_test(OPP)

    deallocate(OPP)
    allocate (t_optprop_3_6 :: OPP)
    call this_test(OPP)
    deallocate(OPP)

    contains
      subroutine this_test(OPP)
        class(t_optprop)  :: OPP

    PETSC_COMM_WORLD = comm
    call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)
    call init_mpi_data_parameters(comm)

    call read_commandline_options()

    call OPP%init([zero],[60._ireals], comm)

    allocate(dir2dir(OPP%OPP_LUT%dir_streams**2))
    allocate(dir2diff(OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams))

    do i=1,ubound(dir2dir,1)
      dir2dir(i) = i
    enddo

    do i=1,ubound(dir2diff,1)
      dir2diff(i) = i
    enddo


    call OPP%dir2dir_coeff_symmetry(dir2dir, .True., .True.)
    call OPP%dir2dir_coeff_symmetry(dir2dir, .True., .True.)

    do i=1,ubound(dir2dir,1)
      @assertEqual(i,dir2dir(i), 'Coeff dir2dir not equal after switching two time north-south and east-west')
    end do

    call OPP%dir2diff_coeff_symmetry(dir2diff, .True., .True.)
    call OPP%dir2diff_coeff_symmetry(dir2diff, .True., .True.)

    do i=1,ubound(dir2diff,1)
      @assertEqual(i,dir2diff(i), 'Coeff dir2diff not equal after switching two time north-south and east-west')
    end do

    deallocate(dir2dir)
    deallocate(dir2diff)
    call OPP%destroy()

    call PetscFinalize(ierr)
    end subroutine
  end subroutine


  @test(npes =[1])
  subroutine test_pprts_symmetry_ex1(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: numnodes, comm
    integer(iintegers) :: myid

    integer(iintegers),parameter :: nxp=9,nyp=9,nv=100
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=10, theta0=60
    real(ireals),parameter :: albedo=0., dz=dx
    real(ireals),parameter :: incSolar=1000
    real(ireals),parameter :: atolerance=1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir0,fdn0,fup0,fdiv0

    real(ireals),allocatable,dimension(:,:,:) :: fdir1,fdn1,fup1,fdiv1

    integer(iintegers) :: i,j,k, ni,nj

    !type(t_solver_3_6)  :: solver
    type(t_solver_8_10) :: solver

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    !!!!!!!!!!!!!!!!!!!!  calculation for phi0 = 0  !!!!!!!!!!!!!!!!!!!!!!!
    call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    allocate(fdir0 (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdn0  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fup0  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdiv0 (solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    allocate(fdir1 (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdn1  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fup1  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdiv1 (solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))


    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = 1._ireals/nv/dz
    ksca = 1._ireals/nv/dz
    g    = zero

    kabs(nv/2,nxp/2+1,nyp/2+1) = 1/dz
    ksca(nv/2,nxp/2+1,nyp/2+1) = 1/dz
    g   (nv/2,nxp/2+1,nyp/2+1) = .9

    call set_optical_properties(solver, albedo, kabs, ksca, g)!, B )
    call set_angles(solver, 10._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=10)

    call set_angles(solver, 80._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=80)

    call set_angles(solver, 100._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=100)

    call set_angles(solver, 190._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=190)

    call set_angles(solver, 280._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=280)

    call pprts_get_result(solver, fdir0,fdn0,fup0,fdiv0, opt_solution_uid=10)
    call pprts_get_result(solver, fdir1,fdn1,fup1,fdiv1, opt_solution_uid=190)

    do j=lbound(fdir0,3), ubound(fdir0,3)
      do i=lbound(fdir0,2), ubound(fdir0,2)
        ni = ubound(fdir0,2)-i+lbound(fdir0,2)
        nj = ubound(fdir0,3)-j+lbound(fdir0,3)

        do k=lbound(fdiv0,1), ubound(fdiv0,1)
          @assertEqual(fdiv0(k,ni,nj), fdiv1(k,i,j), atolerance, '10 -> 190: divergence not symmetric for azimuth')
        enddo
        do k=lbound(fdir0,1), ubound(fdir0,1)
          @assertEqual(fdir0(k,ni,nj), fdir1(k,i,j), atolerance, '10 -> 190: Edirradiation not symmetric for azimuth')
          @assertEqual(fdn0 (k,ni,nj), fdn1 (k,i,j), atolerance, '10 -> 190: Edn radiation not symmetric for azimuth')
          @assertEqual(fup0 (k,ni,nj), fup1 (k,i,j), atolerance, '10 -> 190: Eup radiation not symmetric for azimuth')
        enddo
      enddo
    enddo

    call pprts_get_result(solver, fdir1,fdn1,fup1,fdiv1, opt_solution_uid=80)

    do j=lbound(fdir0,3), ubound(fdir0,3)
      do i=lbound(fdir0,2), ubound(fdir0,2)
        ni = j
        nj = i

        do k=lbound(fdiv0,1), ubound(fdiv0,1)
          @assertEqual(fdiv0(k,ni,nj), fdiv1(k,i,j), atolerance, '10 -> 80: divergence not symmetric for azimuth')
        enddo
        do k=lbound(fdir0,1), ubound(fdir0,1)
          @assertEqual(fdir0(k,ni,nj), fdir1(k,i,j), atolerance, '10 -> 80: Edirradiation not symmetric for azimuth')
          @assertEqual(fdn0 (k,ni,nj), fdn1 (k,i,j), atolerance, '10 -> 80: Edn radiation not symmetric for azimuth')
          @assertEqual(fup0 (k,ni,nj), fup1 (k,i,j), atolerance, '10 -> 80: Eup radiation not symmetric for azimuth')
        enddo
      enddo
    enddo

    call pprts_get_result(solver, fdir0,fdn0,fup0,fdiv0, opt_solution_uid=100)

    !do j=lbound(fdir0,3), ubound(fdir0,3)
    !  do i=lbound(fdir0,2), ubound(fdir0,2)
    !    ni = j
    !    nj = i
    !    nj = ubound(fdir0,3)-nj+lbound(fdir0,3)

    !    do k=lbound(fdiv0,1), ubound(fdiv0,1)
    !      @assertEqual(fdiv0(k,ni,nj), fdiv1(k,i,j), atolerance, '90: divergence not symmetric for azimuth')
    !    enddo
    !    do k=lbound(fdir0,1), ubound(fdir0,1)
    !      @assertEqual(fdir0(k,ni,nj), fdir1(k,i,j), atolerance, '90: Edirradiation not symmetric for azimuth')
    !      @assertEqual(fdn0 (k,ni,nj), fdn1 (k,i,j), atolerance, '90: Edn radiation not symmetric for azimuth')
    !      @assertEqual(fup0 (k,ni,nj), fup1 (k,i,j), atolerance, '90: Eup radiation not symmetric for azimuth')
    !    enddo
    !  enddo
    !enddo

    call pprts_get_result(solver, fdir1,fdn1,fup1,fdiv1, opt_solution_uid=280)

    do j=lbound(fdir0,3), ubound(fdir0,3)
      do i=lbound(fdir0,2), ubound(fdir0,2)
        ni = ubound(fdir0,2)-i+lbound(fdir0,2)
        nj = ubound(fdir0,3)-j+lbound(fdir0,3)

        do k=lbound(fdiv0,1), ubound(fdiv0,1)
          @assertEqual(fdiv0(k,ni,nj), fdiv1(k,i,j), atolerance, '270: divergence not symmetric for azimuth')
        enddo
        do k=lbound(fdir0,1), ubound(fdir0,1)
          @assertEqual(fdir0(k,ni,nj), fdir1(k,i,j), atolerance, '270: Edirradiation not symmetric for azimuth')
          @assertEqual(fdn0 (k,ni,nj), fdn1 (k,i,j), atolerance, '270: Edn radiation not symmetric for azimuth')
          @assertEqual(fup0 (k,ni,nj), fup1 (k,i,j), atolerance, '270: Eup radiation not symmetric for azimuth')
        enddo
      enddo
    enddo
    call destroy_pprts(solver, .True.)
  end subroutine
end module
