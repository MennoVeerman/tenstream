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

module m_optprop_parameters
      use m_data_parameters,only : ireals,iintegers
      implicit none

      !> \page optprop_parameters Parameters concerning the transport coefficients
      !! You should define some parameters about how and where to find the Look
      !! Up Tables foir the transport coefficients
      !> Have a look at the options in m_optprop_parameters::

      !-----------------------------------------
      !> Define the path to the Lookuptables
      !!
      !!  This has to be a reachable path for rank 0,
      !!  At MIM in Munich please set to
      !!  '/home/opt/cosmo_tica_lib/tenstream/optpropLUT/LUT'
      !!  At ZMAW in Hamburg please set to
      !!  '/scratch/mpi/mpiaes/m300362/tenstream_LUT/LUT'
      !-----------------------------------------

      character(len=300) :: lut_basename='./LUT' !'/scratch/mpi/mpiaes/m300362/tenstream_LUT/LUT'

      !-----------------------------------------
      !> Define the mode to calculate coeffs   -
      !!-----------------------------------------
      !! 0 :: retrieve directly from Lookuptable
      !! 1 :: retrieve from ANN (experimental) ::
      !!        this assumes you actually have a
      !!        ANN that can be queried....
      !!
      integer(iintegers) :: coeff_mode = 0

      !-----------------------------------------
      !- Define the size of the Lookuptables:  -
      !-----------------------------------------
      !
      ! You should not need to change this... but feel free to play around...
      ! interp_mode 1 == nearest neighbour interpolation
      ! interp_mode 2 == linear interpolation

!     integer(iintegers) ,parameter :: Ndz_8_10=4, Nkabs_8_10=2, Nksca_8_10=2, Ng_8_10=3, Nphi_8_10=10, Ntheta_8_10=10, Ndiff_8_10=10, Ndir_8_10=8, interp_mode_8_10=1
!     integer(iintegers) ,parameter :: Ndz_8_10=4, Nkabs_8_10=5, Nksca_8_10=5, Ng_8_10=3, Nphi_8_10=10, Ntheta_8_10=10, Ndiff_8_10=10, Ndir_8_10=8, interp_mode_8_10=1
      integer(iintegers) ,parameter :: Ndz_8_10=20, Nkabs_8_10=20, Nksca_8_10=20, Ng_8_10=3, Nphi_8_10=10, Ntheta_8_10=19, Ndiff_8_10=10, Ndir_8_10=8, interp_mode_8_10=4
!      integer(iintegers) ,parameter :: Ndz_8_10=30, Nkabs_8_10=30, Nksca_8_10=30, Ng_8_10=4, Nphi_8_10=10, Ntheta_8_10=10, Ndiff_8_10=10, Ndir_8_10=8, interp_mode_8_10=2
!      integer(iintegers) ,parameter :: Ndz_8_10=50, Nkabs_8_10=80, Nksca_8_10=80, Ng_8_10=5, Nphi_8_10=10, Ntheta_8_10=10, Ndiff_8_10=10, Ndir_8_10=8, interp_mode_8_10=2

!      integer(iintegers) ,parameter :: Ndz_8_10=3, Nkabs_8_10=10, Nksca_8_10=10, Ng_8_10=3, Nphi_8_10=10, Ntheta_8_10=10, interp_mode_8_10=2
!      integer(iintegers) ,parameter :: Ndz_8_10=3, Nkabs_8_10=100, Nksca_8_10=100, Ng_8_10=3, Nphi_8_10=10, Ntheta_8_10=10, interp_mode_8_10=1
!      integer(iintegers) ,parameter :: Ndz_8_10=3, Nkabs_8_10=30, Nksca_8_10=30, Ng_8_10=3, Nphi_8_10=10, Ntheta_8_10=10, interp_mode_8_10=1
!      integer(iintegers) ,parameter :: Ndz_8_10=3, Nkabs_8_10=200, Nksca_8_10=200, Ng_8_10=1, Nphi_8_10=10, Ntheta_8_10=10, interp_mode_8_10=2
!      integer(iintegers) ,parameter :: Ndz_8_10=3, Nkabs_8_10=40, Nksca_8_10=40, Ng_8_10=10, Nphi_8_10=10, Ntheta_8_10=10, interp_mode_8_10=2

      integer(iintegers) ,parameter :: Ndz_1_2=80, Nkabs_1_2=20, Nksca_1_2=20, Ng_1_2=3, Nphi_1_2=1, Ntheta_1_2=19, Ndiff_1_2=2, Ndir_1_2=1, interp_mode_1_2=2

      integer(iintegers) ,parameter :: Ntau=50, Nw0=10

      !-----------------------------------------
      !- Define precision of coefficients      -
      !-----------------------------------------
      ! absolute tolerance and relatice tolerance have to be reached for every
      ! coefficient

!      real(ireals),parameter :: stddev_atol=1e-2_ireals
!      real(ireals),parameter :: stddev_atol=5e-3_ireals
!      real(ireals),parameter :: stddev_atol=2e-3_ireals
      real(ireals),parameter :: stddev_atol=1e-4_ireals
!      real(ireals),parameter :: stddev_atol=5e-6_ireals

!      real(ireals),parameter :: stddev_rtol=5e-2_ireals
      real(ireals),parameter :: stddev_rtol=1e-3_ireals
!      real(ireals),parameter :: stddev_rtol=1e-4_ireals

      ! Do some sanity checks on coefficients -- only disable if you are sure
      ! what to expect.
      logical,parameter :: ldebug_optprop=.False.
!      logical,parameter :: ldebug_optprop=.True.

      ! Use delta scaling on optical properties? -- this significantly reduces
      ! the size of the lookuptables.
      logical,parameter :: ldelta_scale=.True.

      ! Treat direct2diffuse radiation in a cone around solar angle as direct
      ! radiation.
!      real(ireals),parameter :: delta_scale_truncate=.9848_ireals ! .9848 = 10 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling
      real(ireals),parameter :: delta_scale_truncate=1.0_ireals   !1.     = 0 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.8660_ireals ! .8660 = 30 degrees delta scaling


end module
