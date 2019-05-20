!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_optprop

#ifdef _XLF
      use ieee_arithmetic
#define isnan ieee_is_nan
#endif

use m_optprop_parameters, only : ldebug_optprop, coeff_mode, wedge_sphere_radius, param_eps
use m_helper_functions, only : rmse, CHKERR, itoa, ftoa, approx, deg2rad, rad2deg, swap
use m_data_parameters, only: ireals,irealLUT,irealLUT,iintegers,one,zero,i0,i1,inil,mpiint
use m_optprop_LUT, only : t_optprop_LUT, t_optprop_LUT_1_2,t_optprop_LUT_3_6, t_optprop_LUT_3_10, &
  t_optprop_LUT_8_10, t_optprop_LUT_3_16, t_optprop_LUT_8_16, t_optprop_LUT_8_18, &
  t_optprop_LUT_wedge_5_8, t_optprop_LUT_wedge_18_8
use m_optprop_ANN, only : ANN_init, ANN_get_dir2dir, ANN_get_dir2diff, ANN_get_diff2diff
use m_boxmc_geometry, only : setup_default_unit_cube_geometry, setup_default_wedge_geometry
use m_eddington, only: eddington_coeff_zdun
use m_tenstream_options, only: twostr_ratio

use m_LUT_param_phi, only: theta_from_param_theta, iterative_phi_theta_from_param_phi_and_param_theta

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

implicit none

private
public :: t_optprop, t_optprop_cube, t_optprop_wedge, &
  t_optprop_1_2, t_optprop_3_6, t_optprop_3_10, &
  t_optprop_wedge_5_8, t_optprop_wedge_18_8, &
  t_optprop_8_10, t_optprop_3_16, t_optprop_8_16, t_optprop_8_18, &
  OPP_1D_RETCODE, OPP_TINYASPECT_RETCODE

type,abstract :: t_optprop
  logical :: optprop_debug=ldebug_optprop
  class(t_optprop_LUT), allocatable :: OPP_LUT
  contains
    procedure :: init
    procedure :: get_coeff_bmc
    procedure :: destroy
end type

! Cube types
type,abstract,extends(t_optprop) :: t_optprop_cube
  contains
    procedure :: get_coeff => get_coeff_cube
    procedure :: dir2dir_coeff_symmetry => dir2dir_coeff_symmetry_none
    procedure :: dir2diff_coeff_symmetry => dir2diff_coeff_symmetry_none
end type

! we introduce one special cube type for 8 direct streams, this way, all of them can share dir2dir_coeff_symmetry
type,abstract,extends(t_optprop_cube) :: t_optprop_cube_dir8
  contains
    procedure :: dir2dir_coeff_symmetry => dir2dir8_coeff_symmetry
end type

type,extends(t_optprop_cube) :: t_optprop_1_2
end type

type,extends(t_optprop_cube) :: t_optprop_3_6
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff6_coeff_symmetry
end type

type,extends(t_optprop_cube) :: t_optprop_3_10
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff10_coeff_symmetry
end type

type,extends(t_optprop_cube) :: t_optprop_3_16
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff16_coeff_symmetry
end type

type,extends(t_optprop_cube_dir8) :: t_optprop_8_10
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff10_coeff_symmetry
end type

type,extends(t_optprop_cube_dir8) :: t_optprop_8_16
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff16_coeff_symmetry
end type

type,extends(t_optprop_cube_dir8) :: t_optprop_8_18
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff18_coeff_symmetry
end type

! Wedge types
type,abstract,extends(t_optprop) :: t_optprop_wedge
  contains
    procedure :: get_coeff => get_coeff_wedge
end type

type,extends(t_optprop_wedge) :: t_optprop_wedge_5_8
end type

type,extends(t_optprop_wedge) :: t_optprop_wedge_18_8
end type

integer(mpiint), parameter :: OPP_1D_RETCODE = -1_mpiint
integer(mpiint), parameter :: OPP_TINYASPECT_RETCODE = -2_mpiint

contains

  subroutine init(OPP, comm, skip_load_LUT)
      class(t_optprop), intent(inout) :: OPP
      integer(mpiint) ,intent(in) :: comm
      logical, intent(in), optional :: skip_load_LUT
      integer(mpiint) :: ierr

      select case (coeff_mode)
          case(i0) ! LookUpTable Mode
            select type(OPP)
              class is (t_optprop_1_2)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_1_2::OPP%OPP_LUT)

              class is (t_optprop_3_6)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_6::OPP%OPP_LUT)

              class is (t_optprop_3_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_10::OPP%OPP_LUT)

              class is (t_optprop_8_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_10::OPP%OPP_LUT)

              class is (t_optprop_3_16)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_16::OPP%OPP_LUT)

              class is (t_optprop_8_16)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_16::OPP%OPP_LUT)

              class is (t_optprop_8_18)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_18::OPP%OPP_LUT)

              class is (t_optprop_wedge_5_8)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_wedge_5_8::OPP%OPP_LUT)

              class is (t_optprop_wedge_18_8)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_wedge_18_8::OPP%OPP_LUT)
              class default
                call CHKERR(1_mpiint, ' init optprop : unexpected type for optprop object!')
            end select
            call OPP%OPP_LUT%init(comm, skip_load_LUT)

          case(i1) ! ANN
            call ANN_init(comm, ierr)
  !          stop 'ANN not yet implemented'
          case default
            call CHKERR(1_mpiint, 'coeff mode optprop initialization not defined')
        end select
  end subroutine
  subroutine destroy(OPP)
      class(t_optprop) :: OPP
      if(allocated(OPP%OPP_LUT)) then
          call OPP%OPP_LUT%destroy()
          deallocate(OPP%OPP_LUT)
      endif
  end subroutine

  subroutine get_coeff_wedge(OPP, tauz, w0, g, aspect_zx, ldir, C, ierr, wedge_coords, angles)
    class(t_optprop_wedge)              :: OPP
    logical,intent(in)                  :: ldir
    real(irealLUT),intent(in)           :: tauz, w0, g, aspect_zx
    real(irealLUT),intent(in)           :: wedge_coords(:) ! 2 coordinates of wedge C_point, only used for wedge OPP types
    real(irealLUT),intent(in),optional  :: angles(:)
    real(irealLUT),intent(out)          :: C(:)
    integer(mpiint), intent(out) :: ierr

    logical,parameter :: compute_coeff_online=.False.

    ierr = 0

    if(ldebug_optprop) then
      call check_inp(OPP, tauz, w0, g, aspect_zx, ldir, C, angles)
    endif

    if(ldebug_optprop) then
      if(.not.approx(g,0._irealLUT)) call CHKERR(1_mpiint, 'wedge LUT does not have values for other than g==0')
    endif

    if(handle_aspect_zx_1D_case()) then
      return
    endif

    if(ldir .and. compute_coeff_online) then
      call do_bmc_computation(C)
      return
    endif

    call do_wedge_lookup(tauz, w0, aspect_zx, ldir, angles)

    if(.False. .and. ldir) call print_coeff_diff()

    !if(ldir .and. present(angles)) then
    !  call handle_critical_azimuth()
    !endif
    contains
      subroutine do_bmc_computation(Cbmc)
        real(irealLUT), intent(out) :: Cbmc(:)
        real(ireals), allocatable :: vertices(:)
        real(irealLUT) :: phi, theta

        call setup_default_wedge_geometry(&
          real([0,0], ireals), &
          real([1,0], ireals), &
          real(wedge_coords, ireals), &
          real(aspect_zx, ireals), &
          vertices, &
          real(wedge_sphere_radius, ireals))

        call iterative_phi_theta_from_param_phi_and_param_theta(&
          real(vertices, irealLUT), &
          angles(1), angles(2), phi, theta, ierr); call CHKERR(ierr)

        print *,'Cbmc', tauz, w0, g, aspect_zx, wedge_coords, angles, rad2deg(phi), rad2deg(theta)
        call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), ldir, Cbmc, &
          [rad2deg(phi), rad2deg(theta)])
      end subroutine
      subroutine print_coeff_diff()
        real(irealLUT) :: Cbmc(size(C))
        real(ireals) :: err(2)
        integer(iintegers) :: isrc

        call do_bmc_computation(Cbmc)

        err = rmse(real(C, ireals), real(Cbmc, ireals))
        print *,'rmse', err
        do isrc=1,OPP%OPP_LUT%dir_streams
          print *, 'lut src', isrc, ':', C(isrc:OPP%OPP_LUT%dir_streams**2:OPP%OPP_LUT%dir_streams)
          print *, 'bmc src', isrc, ':', Cbmc(isrc:OPP%OPP_LUT%dir_streams**2:OPP%OPP_LUT%dir_streams)
        enddo
        if(err(2).gt.one) then
          !call CHKERR(1_mpiint, 'DEBUG')
        endif
        !C = Cbmc
      end subroutine
      subroutine do_wedge_lookup(tauz, w0, aspect_zx, ldir, angles)
        real(irealLUT), intent(in) :: tauz, w0, aspect_zx
        logical,intent(in)       :: ldir
        real(irealLUT),intent(in),optional :: angles(:)
        real(irealLUT) :: save_param_phi, save_param_theta

          select case (coeff_mode)
          case(i0) ! LookUpTable Mode
            if(present(angles)) then ! obviously we want the direct coefficients

              associate(&
                  Cx => wedge_coords(1), &
                  Cy => wedge_coords(2), &
                  param_phi => angles(1), &
                  param_theta => angles(2) )

                call handle_critical_param_phi(param_phi, save_param_phi)
                call handle_critical_param_theta(param_theta, save_param_theta)

                if(ldir) then ! dir2dir
                  call OPP%OPP_LUT%LUT_get_dir2dir([tauz, w0, aspect_zx, Cx, Cy, &
                    save_param_phi, save_param_theta], C)
                else ! dir2diff
                  call OPP%OPP_LUT%LUT_get_dir2diff([tauz, w0, aspect_zx, Cx, Cy, &
                    save_param_phi, save_param_theta], C)
                endif
              end associate
            else
              ! diff2diff
              associate(Cx => wedge_coords(1), Cy => wedge_coords(2))
                call OPP%OPP_LUT%LUT_get_diff2diff([tauz, w0, aspect_zx, Cx, Cy], C)
              end associate
            endif

          case default
            call CHKERR(1_mpiint, 'particular value of coeff mode in optprop_parameters is not defined: '//itoa(coeff_mode))
          end select
      end subroutine

      subroutine handle_critical_param_phi(param_phi, save_param_phi)
        real(irealLUT), intent(in) :: param_phi
        real(irealLUT), intent(out) :: save_param_phi
        logical :: lsample_critical

        lsample_critical = .False.
        if(approx(abs(param_phi), 1._irealLUT, param_eps)) lsample_critical = .True.

        if(lsample_critical) then
          if(param_phi.le.-1._irealLUT) then
            save_param_phi = -1._irealLUT-param_eps
          elseif(param_phi.ge.1._irealLUT) then !1.0001
            save_param_phi = 1._irealLUT+param_eps
          elseif(param_phi.lt.0._irealLUT) then !-.999
            save_param_phi = -1._irealLUT+param_eps
          else ! .999
            save_param_phi = 1._irealLUT-param_eps
          endif
        else
          save_param_phi = param_phi
        endif
      end subroutine
      subroutine handle_critical_param_theta(param_theta, save_param_theta)
        real(irealLUT), intent(in) :: param_theta
        real(irealLUT), intent(out) :: save_param_theta
        logical :: lsample_critical

        lsample_critical = .False.
        if(approx(abs(param_theta), 0._irealLUT, param_eps)) lsample_critical = .True.

        if(lsample_critical) then
          if(param_theta.lt.0._irealLUT) then
            save_param_theta = -param_eps
          else ! .0001
            save_param_theta = param_eps
          endif
        else
          save_param_theta = param_theta
        endif
      end subroutine

      logical function handle_aspect_zx_1D_case()
        real(ireals) :: c11,c12,c13,c23,c33
        real(irealLUT) :: restricted_aspect_zx
        real(irealLUT) :: mu

        handle_aspect_zx_1D_case = .False.

        !TODO: here we may have incoming radiation at the sides and we just drop that
        ! this has to be fixed for anisotropic grids

        if(present(angles)) then

          if(aspect_zx.ge.twostr_ratio) then
            C = zero
            mu = cos(theta_from_param_theta(angles(2), 0._irealLUT))

            call eddington_coeff_zdun(&
              real(tauz, ireals), &
              real(w0, ireals), &
              real(g, ireals), &
              real(mu, ireals), &
              c11,c12,c13,c23,c33)

            if(ldir) then
              select type(OPP)
              class is (t_optprop_wedge_5_8)
                ! set the transport coeffs for src top to zero, leave the rest.
                C(5*4+1) = real(c33, irealLUT) ! from top to bot
                !C(22:24) = 1 ! from sides to bot
              class is (t_optprop_wedge_18_8)
                C(18*15+1) = real(c33, irealLUT) ! from top to bot
                C(18*16+2) = real(c33, irealLUT) ! from top to bot
                C(18*17+3) = real(c33, irealLUT) ! from top to bot
              class default
                call CHKERR(1_mpiint, 'wedge handle_aspect_zx_1D_case not implemented for this type of OPP')
              end select

            else
              select type(OPP)
              class is (t_optprop_wedge_5_8)
                C(0*5+1) = real(c13, irealLUT) ! reflection
                C(7*5+1) = real(c23, irealLUT) ! transmission
              class is (t_optprop_wedge_18_8)
                C(0*18+1) = real(c13, irealLUT)
                C(0*18+2) = real(c13, irealLUT)
                C(0*18+3) = real(c13, irealLUT)
                C(7*18+1) = real(c23, irealLUT)
                C(7*18+2) = real(c23, irealLUT)
                C(7*18+3) = real(c23, irealLUT)
              class default
                call CHKERR(1_mpiint, 'wedge handle_aspect_zx_1D_case not implemented for this type of OPP')
              end select
            endif

            handle_aspect_zx_1D_case = .True.
            ierr = OPP_1D_RETCODE

          elseif(aspect_zx.lt.OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)) then
            restricted_aspect_zx = min(max(aspect_zx, OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)), &
              OPP%OPP_LUT%dirconfig%dims(3)%vrange(2))
            call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
            handle_aspect_zx_1D_case = .True.
            ierr = OPP_TINYASPECT_RETCODE
          endif

        else ! diffuse

          !if(aspect_zx.gt.OPP%OPP_LUT%diffconfig%dims(3)%vrange(2)) then
          if(aspect_zx.ge.twostr_ratio) then
            C = zero

            call eddington_coeff_zdun(&
              real(tauz, ireals), &
              real(w0, ireals), &
              real(g, ireals), &
              one, &
              c11,c12,c13,c23,c33)

            ! transmission & reflection towards top plate
            C(1)     = real(c12, irealLUT)
            C(8)     = real(c11, irealLUT)
            ! and bot plate
            C(7*8+1) = real(c11, irealLUT)
            C(8*8)   = real(c12, irealLUT)

            handle_aspect_zx_1D_case = .True.

          elseif(aspect_zx.lt.OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)) then
            restricted_aspect_zx = min(max(aspect_zx, OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)), &
              OPP%OPP_LUT%diffconfig%dims(3)%vrange(2))
            call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
            handle_aspect_zx_1D_case = .True.
            ierr = OPP_TINYASPECT_RETCODE
          endif

        endif
      end function
  end subroutine

  subroutine get_coeff_cube(OPP, tauz, w0, g, aspect_zx, dir, C, ierr, angles, lswitch_east, lswitch_north)
    class(t_optprop_cube)             :: OPP
    logical,intent(in)                :: dir
    real(irealLUT),intent(in)           :: tauz, w0, g, aspect_zx
    real(irealLUT),intent(in),optional  :: angles(:)
    logical,intent(in)                  :: lswitch_east, lswitch_north
    real(irealLUT),intent(out)          :: C(:)
    integer(mpiint), intent(out) :: ierr

    logical,parameter :: compute_coeff_online=.False.
    real(ireals), allocatable :: vertices(:)
    real(irealLUT) :: save_aspect_zx
    ierr = 0

    if(compute_coeff_online) then
      call setup_default_unit_cube_geometry(one, one, real(aspect_zx, ireals), vertices)
      call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), dir, C, angles)
      return
    endif

    if(ldebug_optprop) call check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles)


    select case (coeff_mode)

    case(i0) ! LookUpTable Mode

      if(present(angles)) then ! obviously we want the direct coefficients
        if(aspect_zx.lt.OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)) then
          save_aspect_zx = OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)
          ierr = OPP_TINYASPECT_RETCODE
        else
          save_aspect_zx = aspect_zx
        endif
        if(dir) then ! dir2dir
          call OPP%OPP_LUT%LUT_get_dir2dir([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
          call OPP%dir2dir_coeff_symmetry(C, lswitch_east, lswitch_north)
        else         ! dir2diff
          call OPP%OPP_LUT%LUT_get_dir2diff([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
          call OPP%dir2diff_coeff_symmetry(C, lswitch_east, lswitch_north)
        endif
      else
        ! diff2diff
        if(aspect_zx.lt.OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)) then
          save_aspect_zx = OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)
          ierr = OPP_TINYASPECT_RETCODE
        else
          save_aspect_zx = aspect_zx
        endif
        call OPP%OPP_LUT%LUT_get_diff2diff([tauz, w0, save_aspect_zx, g], C)
      endif


    case(i1) ! ANN

      if(present(angles)) then ! obviously we want the direct coefficients
        if(dir) then ! specifically the dir2dir
          call ANN_get_dir2dir(tauz, w0, g, aspect_zx, angles(1), angles(2), C)
        else ! dir2diff
          call ANN_get_dir2diff(tauz, w0, g, aspect_zx, angles(1), angles(2), C)
        endif
      else
        ! diff2diff
        call ANN_get_diff2diff(tauz, w0, g, aspect_zx, C)
      endif

    case default
      call CHKERR(1_mpiint, 'particular value of coeff mode in optprop_parameters is not defined: '//itoa(coeff_mode))
    end select

  end subroutine

  subroutine get_coeff_bmc(OPP, vertices, tauz, w0, g, dir, C, angles)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: tauz, w0, g, vertices(:)
      logical,intent(in) :: dir
      real(irealLUT),intent(out):: C(:)
      real(irealLUT),intent(in),optional :: angles(2)

      real(irealLUT) :: S_diff(OPP%OPP_LUT%diff_streams),T_dir(OPP%OPP_LUT%dir_streams)
      real(irealLUT) :: S_tol (OPP%OPP_LUT%diff_streams),T_tol(OPP%OPP_LUT%dir_streams)
      integer(iintegers) :: isrc

      real(ireals), parameter :: atol=5e-3_ireals, rtol=5e-2_ireals

      if(present(angles)) then
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, &
              real(vertices, ireals), &
              real(tauz, ireals), &
              real(w0, ireals), &
              real(g, ireals), &
              .True.,   &
              real(angles(1), ireals), &
              real(angles(2), ireals), &
              mpi_comm_self, &
              S_diff, T_dir, S_tol, T_tol, &
              inp_atol=atol, inp_rtol=rtol)
            if(dir) then !dir2dir
              C(isrc:OPP%OPP_LUT%dir_streams**2:OPP%OPP_LUT%dir_streams) = T_dir
            else ! dir2diff
              C(isrc:OPP%OPP_LUT%dir_streams*OPP%OPP_LUT%diff_streams:OPP%OPP_LUT%dir_streams) = S_diff
            endif
          enddo
      else
        ! diff2diff
        do isrc=1,OPP%OPP_LUT%diff_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, &
              real(vertices, ireals), &
              real(tauz, ireals), &
              real(w0, ireals), &
              real(g, ireals), &
              .False.,   &
              zero, zero, &
              mpi_comm_self, &
              S_diff, T_dir, S_tol, T_tol, &
              inp_atol=atol, inp_rtol=rtol)
          C(isrc:OPP%OPP_LUT%diff_streams**2:OPP%OPP_LUT%diff_streams) = S_diff
        enddo
      endif ! angles_present

  end subroutine

  subroutine check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles)
    class(t_optprop) :: OPP
    real(irealLUT),intent(in) :: tauz, w0, g, aspect_zx
    logical,intent(in) :: dir
    real(irealLUT),intent(in):: C(:)
    real(irealLUT),intent(in),optional  :: angles(:)
    if( (any([aspect_zx, tauz, w0, g].lt.zero)) .or. (any(isnan([aspect_zx, tauz, w0, g]))) ) then
      print *,'optprop_lookup_coeff :: corrupt optical properties: bg:: ',[aspect_zx, tauz, w0, g]
      call exit
    endif
    !if(.not.approx(g,zero)) then
    !  call CHKERR(1_mpiint, 'currently the LUT calls do not have a dimensions for assym param g. has to be zero')
    !endif
    if(present(angles)) then
      if(dir .and. size(C).ne. OPP%OPP_LUT%dir_streams**2) then
        print *,'direct called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%dir_streams**2
      endif
      if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams) then
        print *,'dir2diffuse called get_coeff with wrong shaped output array:',size(C),'should be',OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams
      endif
    else
      if(dir .and. size(C).ne. OPP%OPP_LUT%diff_streams) then
        print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams
      endif
      if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams**2) then
        print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams**2
      endif
    endif
  end subroutine

  subroutine dir2diff_coeff_symmetry_none(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    return
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  !for solver_3_6 only the offset is changing in those sides which should be switched
  subroutine dir3_to_diff6_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_6)         :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof = 3
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      coeff(7:9)   = newcoeff([1, 2, 3] + dof*3)
      coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
    endif
    if(lswitch_north) then
      newcoeff = coeff
      coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
      coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  !for solver_3_10 the offset is chaning and the destination order
  subroutine dir3_to_diff10_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_10)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof = 3
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      !coeff( 1: 3) = newcoeff([1, 2, 3]        )
      !coeff( 4: 6) = newcoeff([1, 2, 3] + dof*1)
      coeff( 7: 9) = newcoeff([1, 2, 3] + dof*3)
      coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
      coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
      coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
      ! coeff(19:21) = newcoeff([1, 2, 3] + dof*6)
      ! coeff(22:24) = newcoeff([1, 2, 3] + dof*7)
      ! coeff(25:27) = newcoeff([1, 2, 3] + dof*8)
      ! coeff(28:30) = newcoeff([1, 2, 3] + dof*9)
    endif
    if(lswitch_north) then
      newcoeff = coeff
      !coeff( 1: 3) = newcoeff([1, 2, 3]        )
      !coeff( 4: 6) = newcoeff([1, 2, 3] + dof*1)
      !coeff( 7: 9) = newcoeff([1, 2, 3] + dof*2)
      !coeff(10:12) = newcoeff([1, 2, 3] + dof*3)
      !coeff(13:15) = newcoeff([1, 2, 3] + dof*4)
      !coeff(16:18) = newcoeff([1, 2, 3] + dof*5)
      coeff(19:21) = newcoeff([1, 2, 3] + dof*6)
      coeff(22:24) = newcoeff([1, 2, 3] + dof*8)
      coeff(25:27) = newcoeff([1, 2, 3] + dof*7)
      coeff(28:30) = newcoeff([1, 2, 3] + dof*9)
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir3_to_diff16_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_16)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  !for solver_8_10 the offset is chaning and the destination order
  subroutine dir8_to_diff10_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_10)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter :: dof = 8
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      !coeff(1:8)   = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]        )
      !coeff(9:16)  = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*1 )
      coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*3 )
      coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*2 )
      coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*5 )
      coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*4 )
      !coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*6 )
      !coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*7 )
      !coeff(65:72) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*8 )
      !coeff(73:80) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*9 )
    endif
    if (lswitch_north) then
      newcoeff = coeff
      !coeff(1:8)   = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]        )
      !coeff(9:16)  = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*1 )
      !coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*2 )
      !coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*3 )
      !coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*4 )
      !coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*5 )
      coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*6 )
      coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*8 )
      coeff(65:72) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*7 )
      coeff(73:80) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*9 )
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir8_to_diff16_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_16)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  subroutine dir8_to_diff18_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_18)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  subroutine dir2dir_coeff_symmetry_none(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    return
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  subroutine dir2dir8_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube_dir8)   :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof=8
    real(irealLUT)               :: newcoeff(size(coeff))

    if(lswitch_east) then
      newcoeff = coeff
      coeff(1:8)   = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof  )
      coeff(9:16)  = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]      )
      coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*3)
      coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*2)
      coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*4)
      coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*5)
      coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*6)
      coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*7)
    endif
    if (lswitch_north) then
      newcoeff = coeff
      coeff(1:8)   = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*2)
      coeff(9:16)  = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*3)
      coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]      )
      coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof  )
      coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*4)
      coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*5)
      coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*6)
      coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*7)
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine
end module