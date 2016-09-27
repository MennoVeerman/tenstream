@test(npes =[2]) 
subroutine test_tenstream_ex2(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, one, zero, pi, myid

    use m_tenstream, only : init_tenstream, set_angles, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord,C_diff,C_one,C_one1

    use m_tenstream_options, only: read_commandline_options

#include "petsc/finclude/petscdef.h"
    use petsc 

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: iter, k, numnodes, comm
    !integer(iintegers) :: myid

    integer(iintegers),parameter :: nxp=10,nyp=10,nv=10
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: albedo=0, dz=dx
    real(ireals),parameter :: incSolar = 1000
    real(ireals),parameter :: atolerance = 1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,B
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    real(ireals) :: div_target(nv)
    real(ireals) :: dn_target(nv+1)
    real(ireals) :: up_target(nv+1)

    div_target = [-1.875366,  1.6980000E-03,  8.9445002E-03,  9.2916247E-03,  9.6176248E-03,  9.9483747E-03,  1.0287250E-02,  1.0629000E-02,  1.0525875E-02, -0.1125398]
    dn_target  = [0.0     ,  203.0865,  219.3402,  235.7752,  253.1141,  271.3924,  290.6432,  310.9006,  332.1989,  354.5734,  378.0595]
    up_target  = [203.7094,  219.2596,  235.6834,  253.0133,  271.2819,  290.5224,  310.7687,  332.0553,  354.4171,  377.8448,  390.0775]
    
    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()

    call init_tenstream(comm, nv, nxp, nyp, dx, dy, phi0, theta0, dz1d=dz1d)

    allocate(kabs(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(ksca(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(g   (C_one%zm , C_one%xm,  C_one%ym ))
    allocate(B   (C_one1%zm,C_one1%xm, C_one1%ym ))

    ! First solve for solar radiation

    kabs = 1e-8
    ksca = 1e-6
    g    = .5

    B = one

    allocate(fdn  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fup  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fdiv (C_one%zm, C_one%xm, C_one%ym))


    do iter=1,10
      do k=1,5

        call set_optical_properties(albedo, kabs, ksca, g)
        call solve_tenstream(incSolar, opt_solution_uid=k)

        allocate(fdir (C_diff%zm, C_diff%xm, C_diff%ym))
        call tenstream_get_result(fdir, fdn, fup, fdiv, opt_solution_uid=k)
        deallocate(fdir)

      enddo
      do k=6,10

        call set_optical_properties(albedo, kabs, ksca, g, B)
        call solve_tenstream(-one, opt_solution_uid=k)

        call tenstream_get_result(fdir, fdn, fup, fdiv, opt_solution_uid=k)

      enddo
    enddo

    call destroy_tenstream(.True.)

    if(myid.eq.0) then
        print *,'B:', B(:,1,1)*pi
        print *,'kabs:', kabs(:,1,1)
        print *,'fdn:',  fdn (:,1,1)
        print *,'fup:',  fup (:,1,1)
        print *,'fdiv:', fdiv(:,1,1)
    endif

    @assertEqual(div_target, fdiv(:,1,1), atolerance, 'thermal divergence not correct')
    @assertEqual(dn_target,  fdn (:,1,1), atolerance, 'thermal downw flux not correct')
    @assertEqual(up_target,  fup (:,1,1), atolerance, 'thermal upward fl  not correct')

    ! Check that surface emission is the one that we stick in
    @assertEqual(B(ubound(B,1),1,1)*pi, fup (ubound(fup,1),1,1), atolerance, 'Surface Emission not correct')

end subroutine