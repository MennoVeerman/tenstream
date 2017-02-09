@test(npes =[1])
subroutine test_rrtm_lw(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one

!    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
!        tenstream_get_result, getvecpointer, restorevecpointer, &
!        t_coord

    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast

    use m_tenstr_rrtmg, only : tenstream_rrtmg, destroy_tenstream_rrtmg

#include "petsc/finclude/petscdef.h"
    use petsc 

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    integer(iintegers),parameter :: nxp=9, nyp=3, nzp=10 ! local domain size for each rank
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: albedo_th=0, albedo_sol=.3, phi0=180, theta0=60
    real(ireals),parameter :: atolerance = 1

    real(ireals), dimension(nzp+1,nxp,nyp) :: plev, tlev
    real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr 
    real(ireals), dimension(nzp,nxp,nyp) :: lwc, reliq                                         
    real(ireals),allocatable, dimension(:,:,:) :: edir,edn,eup,abso ! [nlev_merged(-1), nxp, nyp]

    character(len=250),parameter :: atm_filename='afglus_100m.dat'

    integer(iintegers) :: i,j,k, nlev, icld
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    logical,parameter :: ldebug=.True.
    logical :: lthermal, lsolar

    PetscErrorCode :: ierr

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    allocate(nxproc(numnodes), source=nxp)   ! domain decomp in x direction with all ranks
    allocate(nyproc(1), source=nyp)          ! no domain decomposition in y direction 

    call init_mpi_data_parameters(comm)

    do k=1,nzp+1
      plev(k,:,:) = 1000_ireals - (k-one)*500._ireals/(nzp)
      tlev(k,:,:) = 288._ireals - (k-one)*50._ireals/(nzp)
    enddo

    h2ovmr = zero
    o3vmr  = zero
    co2vmr = zero
    ch4vmr = zero
    n2ovmr = zero
    o2vmr  = zero

    lwc = 0
    reliq = 0

    icld = (nzp+1)/2
    lwc  (icld, :,:) = 1e-2
    reliq(icld, :,:) = 10

    tlev (icld  , :,:) = 288
    tlev (icld+1, :,:) = tlev (icld  , :,:)

    tlay = (tlev(1:nzp,:,:) + tlev(2:nzp+1,:,:))/2

!    ! For comparison, compute lw and sw separately
!    if(myid.eq.0 .and. ldebug) print *,'Computing Solar Radiation:'
!    lthermal=.False.; lsolar=.True.
!
!    call tenstream_rrtmg(comm, dx, dy, phi0, theta0, albedo_th, albedo_sol,  &
!      atm_filename, lthermal, lsolar,               &
!      edir,edn,eup,abso,                            & 
!      d_plev=plev, d_tlev=tlev, d_tlay=tlay, d_lwc=lwc, d_reliq=reliq,         &
!      nxproc=nxproc, nyproc=nyproc, opt_time=one*k)
!
!    nlev = ubound(edn,1)
!    if(myid.eq.0) then
!        if(ldebug) then
!            do k=1,nlev
!                print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
!            enddo
!        endif
!
!        @assertEqual(313.47, edir(nlev,1,1), atolerance, 'solar at surface :: direct flux not correct')
!        @assertEqual(143.32, edn (nlev,1,1), atolerance, 'solar at surface :: downw flux not correct')
!        @assertEqual(137.04, eup (nlev,1,1), atolerance, 'solar at surface :: upward fl  not correct')
!        @assertEqual(-1.395E-02, abso(nlev-1,1,1), atolerance, 'solar at surface :: absorption not correct')
!
!        @assertEqual(684.1109, edir(1,1,1), atolerance, 'solar at TOA :: direct flux not correct')
!        @assertEqual(0       , edn (1,1,1), atolerance, 'solar at TOA :: downw flux not correct')
!        @assertEqual(207.18  , eup (1,1,1), atolerance, 'solar at TOA :: upward fl  not correct')
!        @assertEqual(-2.063E-04, abso(1,1,1), atolerance, 'solar at TOA :: absorption not correct')
!
!        @assertEqual(502.23 , edir(nlev-icld  ,1,1), atolerance, 'solar at icloud :: direct flux not correct')
!        @assertEqual(339.22 , edir(nlev-icld+1,1,1), atolerance, 'solar at icloud+1 :: direct flux not correct')
!        @assertEqual(143.68 , edn (nlev-icld+1,1,1), atolerance, 'solar at icloud :: downw flux not correct')
!        @assertEqual(190.29 , eup (nlev-icld  ,1,1), atolerance, 'solar at icloud :: upward fl  not correct')
!        @assertEqual(-0.0242, abso(nlev-icld  ,1,1), atolerance, 'solar at icloud :: absorption not correct')
!    endif


    if(myid.eq.0 .and. ldebug) print *,'Computing Thermal Radiation:'
    lthermal=.True.; lsolar=.False.

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0, albedo_th, albedo_sol,  &
      atm_filename, lthermal, lsolar,               &
      edir,edn,eup,abso,                            & 
      d_plev=plev, d_tlev=tlev, d_tlay=tlay, d_lwc=lwc, d_reliq=reliq,         &
      nxproc=nxproc, nyproc=nyproc, opt_time=one*k)

    nlev = ubound(edn,1)
    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlev
                print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
            enddo
        endif

        @assertEqual(143.32, edn (nlev,1,1), atolerance, 'thermal at surface :: downw flux not correct')
        @assertEqual(137.04, eup (nlev,1,1), atolerance, 'thermal at surface :: upward fl  not correct')
        @assertEqual(-1.395E-02, abso(nlev-1,1,1), atolerance, 'thermal at surface :: absorption not correct')

        @assertEqual(0             , edn (1,1,1), atolerance, 'thermal at TOA :: downw flux not correct')
        @assertEqual(207.18      , eup (1,1,1), atolerance, 'thermal at TOA :: upward fl  not correct')
        @assertEqual(-2.063E-04, abso(1,1,1), atolerance, 'thermal at TOA :: absorption not correct')

        @assertEqual(143.68 , edn (nlev-icld+1,1,1), atolerance, 'thermal at icloud :: downw flux not correct')
        @assertEqual(190.29 , eup (nlev-icld  ,1,1), atolerance, 'thermal at icloud :: upward fl  not correct')
        @assertEqual(-0.0242, abso(nlev-icld  ,1,1), atolerance, 'thermal at icloud :: absorption not correct')
    endif

    call destroy_tenstream_rrtmg()

end subroutine
