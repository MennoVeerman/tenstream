module m_example_uvspec_cld_file
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_pprts_base, only : t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, &
    i0, i1, i2, zero, one, default_str_len

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, hydrostat_dp

  use m_helper_functions, only: CHKERR, itoa, imp_bcast, reverse, spherical_2_cartesian, resize_arr, &
    domain_decompose_2d
  use m_netcdfio, only: ncload, ncwrite, get_global_attribute


  implicit none

contains
  subroutine example_uvspec_cld_file(comm, &
      cldfile, atm_filename, outfile, &
      albedo_th, albedo_sol, &
      lsolar, lthermal, &
      phi0, theta0, &
      Tsrfc, dTdz)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: Tsrfc, dTdz

    real(ireals), dimension(:,:,:), allocatable, target :: lwc, reliq, tlay ! will have global shape Nz, Nx, Ny
    real(ireals), dimension(:,:,:), allocatable, target :: lwc_local, reliq_local, plev_local, tlay_local ! will have global shape Nz, Nx, Ny
    real(ireals), dimension(:,:,:), allocatable, target :: plev, tlev ! will have local shape nzp+1, nxp, nyp
    real(ireals), dimension(:,:), allocatable :: tsfc, tsfc_local ! nxp,nyp
    real(ireals), dimension(:), allocatable :: hhl,xx,yy ! im Nz+1
    character(len=default_str_len) :: groups(2)

    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    real(ireals),allocatable, dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf

    class(t_solver), allocatable :: pprts_solver

    real(ireals) :: dx, dy
    integer(mpiint) :: myid, N_ranks_x, N_ranks_y, ierr
    integer(iintegers) :: iproc, jproc, is, ie, js, je
    integer(iintegers) :: nxp, nyp, nzp ! local sizes of domain, nzp being number of layers

    integer(iintegers) :: i,j,k

    call mpi_comm_rank(comm, myid, ierr)

    ! Load LibRadtran Cloud File
    !call get_global_attribute(cldfile, 'dx', dx)
    !call get_global_attribute(cldfile, 'dy', dy)
  if (myid.eq.0) then
    groups(1) = trim(cldfile)
    groups(2) = trim('lwc'); call ncload(groups, lwc, ierr); call CHKERR(ierr)
    groups(2) = trim('reff'); call ncload(groups, reliq, ierr); call CHKERR(ierr)
    groups(2) = trim('plev'); call ncload(groups, plev, ierr); call CHKERR(ierr)
    groups(2) = trim('tabs'); call ncload(groups, tlay, ierr); call CHKERR(ierr)
    groups(2) = trim('tsfc'); call ncload(groups, tsfc, ierr); call CHKERR(ierr)
    groups(2) = trim('z_lev'); call ncload(groups, hhl, ierr); call CHKERR(ierr)
    groups(2) = trim('x'); call ncload(groups, xx, ierr); call CHKERR(ierr)
    groups(2) = trim('y'); call ncload(groups, yy, ierr); call CHKERR(ierr)
    endif

    call imp_bcast(comm,lwc,0_mpiint)
    call imp_bcast(comm,reliq,0_mpiint)
    call imp_bcast(comm,hhl,0_mpiint)
    call imp_bcast(comm,plev,0_mpiint)
    call imp_bcast(comm,tlay,0_mpiint)
    call imp_bcast(comm,tsfc,0_mpiint)
    call imp_bcast(comm,xx,0_mpiint)
    call imp_bcast(comm,yy,0_mpiint)
    if(size(lwc  ,dim=2).eq.1) call resize_arr(3_iintegers, lwc  , dim=2, lrepeat=.True.)
    if(size(reliq,dim=2).eq.1) call resize_arr(3_iintegers, reliq, dim=2, lrepeat=.True.)

    if(size(lwc  ,dim=3).eq.1) call resize_arr(3_iintegers, lwc  , dim=3, lrepeat=.True.)
    if(size(reliq,dim=3).eq.1) call resize_arr(3_iintegers, reliq, dim=3, lrepeat=.True.)

    dx  = xx(2)-xx(1)  !* 1e+3_ireals
    dy  = yy(2)-yy(1)  !* 1e+3_ireals
    hhl = hhl !# 1e+3_ireals

    if(myid.eq.0) then
      print *,'Loaded LibRadtran Cloud File with:'
      print *,'dx, dy:', dx, dy
      print *,'hhl', hhl
      print *,'shape lwc ', shape(lwc)
      print *,'shape reliq', shape(reliq)
    endif

    ! Determine Domain Decomposition
    call domain_decompose_2d(comm, N_ranks_x, N_ranks_y, ierr); call CHKERR(ierr)
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y

    nxp = size(lwc, dim=2) / N_ranks_x
    nyp = size(lwc, dim=3) / N_ranks_y
    call CHKERR(modulo(size(lwc, dim=2,kind=mpiint), int(N_ranks_x,mpiint)), &
      'x-dimension is not evenly distributable on given communicator!'// &
      'cant put'//itoa(size(lwc, dim=2))//' pixels on '//itoa(N_ranks_x)//' ranks')
    call CHKERR(modulo(size(lwc, dim=3, kind=mpiint), int(N_ranks_y, mpiint)), &
      'y-dimension is not evenly distributable on given communicator!'// &
      'cant put'//itoa(size(lwc, dim=3))//' pixels on '//itoa(N_ranks_y)//' ranks')

    if(myid.eq.0) then
      print *,'Local Domain sizes are:',nxp,nyp
    endif
    jproc = myid / N_ranks_x
    iproc = modulo(myid, int(N_ranks_x, mpiint))

    js = 1 + jproc * nyp
    je = js + nyp -1

    is = 1 + (myid - jproc*N_ranks_x) * nxp
    ie = is + nxp -1

    print *,myid,'i,j proc', iproc, jproc,' local portion: ', is, ie, 'and', js, je

    nzp = size(hhl)-1
    allocate(lwc_local(nzp,nxp,nyp))
    allocate(reliq_local(nzp,nxp,nyp))
    allocate(plev_local(nzp+1,nxp,nyp))
    allocate(tlay_local(nzp,nxp,nyp))
    allocate(tsfc_local(nxp,nyp))

    lwc_local = lwc(:, is:ie, js:je)
    reliq_local = reliq(:, is:ie, js:je)
    plev_local = plev(:, is:ie, js:je)
    tlay_local = tlay(:, is:ie, js:je)
    tsfc_local = tsfc(is:ie, js:je)
    
    deallocate(lwc)
    deallocate(reliq)
    deallocate(plev)
    deallocate(tlay)
    deallocate(tsfc)
    call mpi_barrier(comm, ierr)
    reliq_local = min(max(reliq_local(:,:,:), 2.5), 60.)
    
    allocate(tlev(nzp+1, nxp, nyp))
    tlev(1,:,:) = Tsfc_local
    do k=2,nzp
      tlev(k,:,:) = (tlay_local(k-1,:,:) + tlay_local(k,:,:) )/2.
    enddo
    tlev(nzp+1,:,:) = 2.*tlay_local(nzp,:,:) - tlev(nzp,:,:)
    
    if(myid.eq.0) then
      do k=1,nzp+1
        print *, k, 'plev', plev_local(k,1,1), 'Tlev', tlev(k,1,1)
      enddo
    endif

    call allocate_pprts_solver_from_commandline(pprts_solver, default_solver='3_10')

    call run_rrtmg_lw_sw(pprts_solver, atm_filename, dx, dy, phi0, theta0, &
      plev_local, tlev, &
      lwc_local, reliq_local, &
      albedo_th, albedo_sol, lsolar, lthermal, &
      edir, edn, eup, abso)

    groups(1) = trim(outfile)

    if(allocated(edir)) then
      call gather_all_toZero(pprts_solver%C_one_atm1, edir, gedir)
      if(myid.eq.0) then
        print *,'dumping direct radiation with local and global shape', shape(edir), ':', shape(gedir)
        groups(2) = 'edir'; call ncwrite(groups, gedir, ierr); call CHKERR(ierr)
      endif
    endif
    call gather_all_toZero(pprts_solver%C_one_atm1, edn, gedn)
    call gather_all_toZero(pprts_solver%C_one_atm1, eup, geup)
    call gather_all_toZero(pprts_solver%C_one_atm, abso, gabso)
    if(myid.eq.0) then
      print *,'dumping edn radiation with local and global shape', shape(edn), ':', shape(gedn)
      groups(2) = 'edn' ; call ncwrite(groups, gedn , ierr); call CHKERR(ierr)
      print *,'dumping eup radiation with local and global shape', shape(eup), ':', shape(geup)
      groups(2) = 'eup' ; call ncwrite(groups, geup , ierr); call CHKERR(ierr)
      print *,'dumping abso radiation with local and global shape', shape(abso), ':', shape(gabso)
      groups(2) = 'abso'; call ncwrite(groups, gabso, ierr); call CHKERR(ierr)
    endif

    call destroy_pprts_rrtmg(pprts_solver, lfinalizepetsc=.True.)
  end subroutine

  subroutine run_rrtmg_lw_sw(pprts_solver, atm_filename, dx, dy, phi0, theta0, &
      plev, tlev, lwc, reliq, albedo_th, albedo_sol, lsolar, lthermal, &
      edir, edn, eup, abso)
    class(t_solver) :: pprts_solver
    real(ireals),intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in), dimension(:,:,:), contiguous, target :: lwc, reliq ! dim(Nz,Nx,Ny)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lsolar, lthermal ! switches if solar or thermal computations should be done
    real(ireals), dimension(:,:,:), contiguous, target, intent(in) :: plev ! pressure on layer interfaces [hPa]   dim=nzp+1,nxp,nyp
    real(ireals), dimension(:,:,:), contiguous, target, intent(in) :: tlev ! Temperature on layer interfaces [K]  dim=nzp+1,nxp,nyp

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: comm, myid, N_ranks_x, N_ranks_y, ierr


    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    ! real(ireals), dimension(nzp,nxp,nyp), target :: lwc, reliq

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=*), intent(in) :: atm_filename

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq

    logical,parameter :: ldebug=.True.

    type(t_tenstr_atm) :: atm

    comm = mpi_comm_world
    call mpi_comm_rank(comm, myid, ierr)

    ! Determine Domain Decomposition
    call domain_decompose_2d(comm, N_ranks_x, N_ranks_y, ierr); call CHKERR(ierr)
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y

    allocate(nxproc(N_ranks_x), source=size(plev,2, kind=iintegers)) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=size(plev,3, kind=iintegers)) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)
    ! h2ovmr = zero
    ! o3vmr  = zero
    ! co2vmr = zero
    ! ch4vmr = zero
    ! n2ovmr = zero
    ! o2vmr  = zero

    if(myid.eq.0 .and. ldebug) print *,'Setup Atmosphere...'

    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(tlev,1),1:size(tlev,2)*size(tlev,3)) => tlev
    plwc (1:size(lwc ,1),1:size(lwc ,2)*size(lwc ,3)) => lwc
    preliq(1:size(reliq,1),1:size(reliq,2)*size(reliq,3)) => reliq

    call setup_tenstr_atm(comm, .False., atm_filename, &
      pplev, ptlev, atm, &
      d_lwc=plwc, d_reliq=preliq)

    call pprts_rrtmg(comm, pprts_solver, atm, &
      size(plev,2, kind=iintegers), size(plev,3, kind=iintegers), &
      dx, dy, phi0, theta0,                    &
      albedo_th, albedo_sol,                   &
      lthermal, lsolar,                        &
      edir, edn, eup, abso,                    &
      nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn,1)
    if(myid.eq.0) then
      if(ldebug) then
        do k=1,nlev
          if(allocated(edir)) then
          print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
        else
          print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
        endif
        enddo
      endif

      if(allocated(edir)) &
        print *,'surface :: direct flux', edir(nlev,1,1)
      print *,'surface :: downw flux ', edn (nlev,1,1)
      print *,'surface :: upward fl  ', eup (nlev,1,1)
      print *,'surface :: absorption ', abso(nlev-1,1,1)

      if(allocated(edir)) &
        print *,'TOA :: direct flux', edir(1,1,1)
      print *,'TOA :: downw flux ', edn (1,1,1)
      print *,'TOA :: upward fl  ', eup (1,1,1)
      print *,'TOA :: absorption ', abso(1,1,1)

    endif

    ! Tidy up
    call destroy_tenstr_atm(atm)
  end subroutine


end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only : mpiint
  use m_example_uvspec_cld_file

  implicit none

  integer(mpiint) :: ierr, myid
  logical :: lthermal, lsolar, lflg
  character(len=10*default_str_len) :: cldfile, outfile
  real(ireals) :: Ag, phi0, theta0, Tsrfc, dTdz
  character(len=default_str_len) :: atm_filename

  logical :: luse_plexrt
  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-cld', cldfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply a cloud filename... please call with -cld <libRadtran_cloud_file.nc>')

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply a output filename... please call with -out <output.nc>')

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-atm_filename', atm_filename, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply an atmosphere filename... please call with -atm_filename <atm.dat>')

  Ag = .1
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

  lsolar = .True.
  lthermal = .False.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)


  phi0 = 180
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg,ierr) ; call CHKERR(ierr)
  theta0 = 60
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg,ierr) ; call CHKERR(ierr)

!  Tsrfc = 288
!  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Tsrfc", Tsrfc, lflg,ierr) ; call CHKERR(ierr)
  

    call example_uvspec_cld_file(mpi_comm_world, cldfile, atm_filename, outfile, &
      zero, Ag, lsolar, lthermal, phi0, theta0, Tsrfc, dTdz)
  call PetscFinalize(ierr)
  call mpi_finalize(ierr)
end program
