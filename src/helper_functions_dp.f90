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

module m_helper_functions_dp
      use m_data_parameters,only : iintegers,ireal_dp,imp_real_dp,imp_iinteger,imp_logical,mpiint, pi_dp
      use m_helper_functions, only: CHKERR
      use mpi

      implicit none

      private
      public imp_bcast, norm, cross_2d, cross_3d, deg2rad, rad2deg, rmse, mean, approx, rel_approx,&
          delta_scale_optprop, delta_scale, cumsum, inc, swap, &
          mpi_logical_and, mpi_logical_or, imp_allreduce_min, imp_allreduce_max, imp_reduce_sum, &
          pnt_in_triangle, pnt_in_rectangle, compute_normal_3d, hit_plane, spherical_2_cartesian, &
          rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, determine_normal_direction, &
          distance_to_edge, distances_to_triangle_edges, triangle_intersection, square_intersection

      interface imp_bcast
        module procedure imp_bcast_real_1d,imp_bcast_real_2d,imp_bcast_real_3d,imp_bcast_real_5d,imp_bcast_int_1d,imp_bcast_int_2d,imp_bcast_int,imp_bcast_real,imp_bcast_logical
      end interface
      interface swap
        module procedure swap_iintegers, swap_ireal_dp
      end interface

      integer(mpiint) :: mpierr
      real(ireal_dp),parameter :: zero=0, one=1

    contains
      pure elemental subroutine inc(x,i)
          real(ireal_dp),intent(inout) :: x
          real(ireal_dp),intent(in) :: i
          x=x+i
      end subroutine
      pure elemental subroutine swap_iintegers(x,y)
        integer(iintegers),intent(inout) :: x,y
        integer(iintegers) :: tmp
        tmp = x
        x = y
        y = tmp
      end subroutine
      pure elemental subroutine swap_ireal_dp(x,y)
        real(ireal_dp),intent(inout) :: x,y
        real(ireal_dp) :: tmp
        tmp = x
        x = y
        y = tmp
      end subroutine

      pure function norm(v)
        real(ireal_dp) :: norm
        real(ireal_dp),intent(in) :: v(:)
        norm = sqrt(dot_product(v,v))
      end function

      pure function cross_3d(a, b)
        real(ireal_dp), dimension(3), intent(in) :: a, b
        real(ireal_dp), dimension(3) :: cross_3d

        cross_3d(1) = a(2) * b(3) - a(3) * b(2)
        cross_3d(2) = a(3) * b(1) - a(1) * b(3)
        cross_3d(3) = a(1) * b(2) - a(2) * b(1)
      end function cross_3d

      pure function cross_2d(a, b)
        real(ireal_dp), dimension(2), intent(in) :: a, b
        real(ireal_dp) :: cross_2d

        cross_2d = a(1) * b(2) - a(2) * b(1)
      end function cross_2d

      elemental function deg2rad(deg)
          real(ireal_dp) :: deg2rad
          real(ireal_dp),intent(in) :: deg
          deg2rad = deg * pi_dp / 180
      end function
      elemental function rad2deg(rad)
        real(ireal_dp) :: rad2deg
        real(ireal_dp),intent(in) :: rad
        rad2deg = rad / pi_dp * 180
      end function

      pure function rmse(a,b)
          real(ireal_dp) :: rmse(2)
          real(ireal_dp),intent(in) :: a(:),b(:)
          rmse(1) = sqrt( mean( (a-b)**2 ) )
          rmse(2) = rmse(1)/max( mean(b), epsilon(rmse) )
      end function

      pure function mean(arr)
          real(ireal_dp) :: mean
          real(ireal_dp),intent(in) :: arr(:)
          mean = sum(arr)/size(arr)
      end function

      elemental logical function approx(a,b,precision)
          real(ireal_dp),intent(in) :: a,b
          real(ireal_dp),intent(in),optional :: precision
          real(ireal_dp) :: factor
          if(present(precision) ) then
            factor = precision
          else
            factor = 10 * epsilon(b)
          endif
          if( a.le.b+factor .and. a.ge.b-factor ) then
            approx = .True.
          else
            approx = .False.
          endif
      end function
      elemental logical function rel_approx(a,b,precision)
          real(ireal_dp),intent(in) :: a,b
          real(ireal_dp),intent(in),optional :: precision
          real(ireal_dp) :: factor,rel_error
          if(present(precision) ) then
            factor = precision
          else
            factor = 10*epsilon(b)
          endif
          rel_error = abs( (a-b)/ max(epsilon(a), ( (a+b)*.5_ireal_dp ) ) )

          if( rel_error .lt. precision ) then
            rel_approx = .True.
          else
            rel_approx = .False.
          endif
      end function


      function mpi_logical_and(comm, lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_and
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_and, 1_mpiint, imp_logical, MPI_LAND, comm, mpierr); call CHKERR(mpierr)
      end function
      function mpi_logical_or(comm, lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_or
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_or, 1_mpiint, imp_logical, MPI_LOR, comm, mpierr); call CHKERR(mpierr)
      end function

      subroutine imp_allreduce_min(comm, v,r)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(in) :: v
          real(ireal_dp),intent(out) :: r
          call mpi_allreduce(v,r,1,imp_real_dp, MPI_MIN,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_allreduce_max(comm, v,r)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(in) :: v
          real(ireal_dp),intent(out) :: r
          call mpi_allreduce(v,r,1,imp_real_dp, MPI_MAX,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_reduce_sum(comm, v, myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(inout) :: v
          integer(mpiint),intent(in) :: myid
          integer(mpiint) :: commsize

          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return 

          if(myid.eq.0) then
            call mpi_reduce(MPI_IN_PLACE, v, 1, imp_real_dp, MPI_SUM, 0, comm, mpierr); call CHKERR(mpierr)
          else
            call mpi_reduce(v, MPI_IN_PLACE, 1, imp_real_dp, MPI_SUM, 0, comm, mpierr); call CHKERR(mpierr)
          endif
      end subroutine

      subroutine  imp_bcast_logical(comm, val, sendid)
          integer(mpiint),intent(in) :: comm
          logical,intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val, 1_mpiint, imp_logical, sendid, comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int(comm, val,sendid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val,1_mpiint,imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_1d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot

          if(sendid.eq.myid) Ntot = size(arr)
          call mpi_bcast(Ntot,1_mpiint,imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_2d(comm, arr,sendid,myid)!
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(2)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real(comm, val,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid

          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val,1_mpiint,imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_1d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = size(arr)
          call imp_bcast_int(comm, Ntot, sendid)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_2d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(2)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_3d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(3)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,3_mpiint,imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3) ) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_5d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:,:,:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(5)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,5_mpiint,imp_iinteger,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3), Ntot(4), Ntot(5) ) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine

      elemental subroutine delta_scale( kabs,ksca,g,factor )
          real(ireal_dp),intent(inout) :: kabs,ksca,g ! kabs, ksca, g
          real(ireal_dp),intent(in),optional :: factor
          real(ireal_dp) :: dtau, w0
          dtau = max( kabs+ksca, epsilon(dtau) )
          w0   = ksca/dtau
          g    = g

          if(present(factor)) then
            call delta_scale_optprop( dtau, w0, g, factor)
          else
            call delta_scale_optprop( dtau, w0, g)
          endif

          kabs= dtau * (one-w0)
          ksca= dtau * w0
      end subroutine
      elemental subroutine delta_scale_optprop( dtau, w0, g,factor)
          real(ireal_dp),intent(inout) :: dtau,w0,g
          real(ireal_dp),intent(in),optional :: factor
          real(ireal_dp) :: f

          g = min( g, one-epsilon(g)*10)
          if(present(factor)) then
            f = factor
          else
            f = g**2
          endif
          dtau = dtau * ( one - w0 * f )
          g    = ( g - f ) / ( one - f )
          w0   = w0 * ( one - f ) / ( one - f * w0 )
      end subroutine

      function cumsum(arr)
          real(ireal_dp),intent(in) :: arr(:)
          real(ireal_dp) :: cumsum(size(arr))
          integer :: i
          cumsum(1) = arr(1)
          do i=2,size(arr)
            cumsum(i) = cumsum(i-1) + arr(i)
          enddo
      end function

      pure function compute_normal_3d(p1,p2,p3)
        ! for a triangle p1, p2, p3, if the vector U = p2 - p1 and the vector V = p3 - p1 then the normal
        ! N = U X V and can be calculated by:
        real(ireal_dp), intent(in) :: p1(3), p2(3), p3(3)
        real(ireal_dp) :: compute_normal_3d(3)
        real(ireal_dp) :: U(3), V(3)

        U = p2-p1
        V = p3-p1

        compute_normal_3d(1) = U(2)*V(3) - U(3)*V(2)
        compute_normal_3d(2) = U(3)*V(1) - U(1)*V(3)
        compute_normal_3d(3) = U(1)*V(2) - U(2)*V(1)

        compute_normal_3d = compute_normal_3d / norm(compute_normal_3d)
      end function

    !> @brief For local azimuth and zenith angles, return the local cartesian vectors phi azimuth, theta zenith angles, angles are input in degrees.
    !> @details theta == 0 :: z = -1, i.e. downward
    !> @details azimuth == 0 :: vector going toward minus y, i.e. sun shines from the north
    !> @details azimuth == 90 :: vector going toward minus x, i.e. sun shines from the east
    pure function spherical_2_cartesian(phi, theta, r)
      real(ireal_dp), intent(in) :: phi, theta
      real(ireal_dp), intent(in), optional :: r

      real(ireal_dp) :: spherical_2_cartesian(3)

      spherical_2_cartesian(1) = -sin(deg2rad(theta)) * sin(deg2rad(phi))
      spherical_2_cartesian(2) = -sin(deg2rad(theta)) * cos(deg2rad(phi))
      spherical_2_cartesian(3) = -cos(deg2rad(theta))

      if(present(r)) spherical_2_cartesian = spherical_2_cartesian*r
    end function

    pure function angle_between_two_vec(p1,p2)
      real(ireal_dp),intent(in) :: p1(:), p2(:)
      real(ireal_dp) :: angle_between_two_vec
      real(ireal_dp) :: n1, n2
      n1 = norm(p1)
      n2 = norm(p2)
      angle_between_two_vec = acos(dot_product(p1/norm(p1), p2/norm(p2)))
    end function

      !> @brief determine distance where a photon p intersects with a plane
      !> @details inputs are the location and direction of a photon aswell as the origin and surface normal of the plane
      pure function hit_plane(p_loc, p_dir, po, pn)
        real(ireal_dp) :: hit_plane
        real(ireal_dp),intent(in) :: p_loc(3), p_dir(3)
        real(ireal_dp),intent(in) :: po(3), pn(3)
        real(ireal_dp) :: discr
        discr = dot_product(p_dir,pn)
        if( ( discr.le. epsilon(discr) ) .and. ( discr.gt.-epsilon(discr)  ) ) then
          hit_plane = huge(hit_plane)
        else
          hit_plane = dot_product(po-p_loc, pn) / discr
        endif
      end function

    subroutine square_intersection(origin, direction, tA, tB, tC, tD, lhit, hit)
      real(ireal_dp), intent(in) :: origin(:), direction(:), tA(:), tB(:), tC(:), tD(:)
      logical, intent(out) :: lhit
      real(ireal_dp), intent(out) :: hit(:)

      lhit = .False.
      hit = huge(hit)

      ! 2 Triangles incorporating, cut along (AC)
      call triangle_intersection(origin, direction, tA, tC, tB, lhit, hit)
      if(lhit) return
      call triangle_intersection(origin, direction, tA, tC, tD, lhit, hit)
      if(lhit) return
    end subroutine

    ! Watertight ray -> triangle intersection code from http://jcgt.org/published/0002/01/05/
    subroutine triangle_intersection(origin, direction, tA, tB, tC, lhit, hit)
      real(ireal_dp), intent(in) :: origin(:), direction(:), tA(:), tB(:), tC(:)
      logical, intent(out) :: lhit
      real(ireal_dp), intent(out) :: hit(:)

      logical, parameter :: ldebug = .False. , BACKFACE_CULLING=.False., HIT_EDGE=.True.

      real(ireal_dp) :: org(0:2), dir(0:2), A(0:2), B(0:2), C(0:2)
      integer(iintegers) :: kx, ky, kz
      real(ireal_dp) :: Sx, Sy, Sz
      real(ireal_dp) :: Ax, Ay, Bx, By, Cx, Cy
      real(ireal_dp) :: Az, Bz, Cz, T
      real(ireal_dp) :: U, V, W
      real(ireal_dp) :: b0, b1, b2
      real(ireal_dp) :: det, rcpDet, dist_times_det


      real(ireal_dp) :: CxBy, CyBx, AxCy, AyCx, BxAy, ByAx

      lhit = .True.
      hit = huge(hit)

      org = origin
      dir = direction

      if(ldebug) print *,'initial direction:', dir
      if(ldebug) print *,'initial origin   :', origin
      if(ldebug) print *,'Triangle coord   :', tA
      if(ldebug) print *,'Triangle coord   :', tB
      if(ldebug) print *,'Triangle coord   :', tC
      ! calculate dimension where the ray direction is maximal (C indexing)
      kz = maxloc(abs(dir), dim=1)-1
      kx = kz+1; if (kx == 3) kx = 0
      ky = kx+1; if (ky == 3) ky = 0
      if(ldebug) print *,'max direction:', kx, ky, kz

      ! swap kx and ky dimension to preserve winding direction of triangles
      if (dir(kz) < zero) call swap(kx, ky)
      if(ldebug) print *,'max direction after swap:', kx, ky, kz
      if(ldebug) print *,'principal direction:', dir(kx), dir(ky), dir(kz)

      ! calculate shear constants
      Sx = dir(kx) / dir(kz)
      Sy = dir(ky) / dir(kz)
      Sz = one / dir(kz)
      if(ldebug) print *,'Shear constants:', Sx, Sy, Sz

      ! calculate vertices relative to ray origin
      A = tA-origin
      B = tB-origin
      C = tC-origin
      if(ldebug) print *,'relative Triangle coords A:', A
      if(ldebug) print *,'relative Triangle coords B:', B
      if(ldebug) print *,'relative Triangle coords C:', C


      ! perform shear and scale of vertices
      Ax = A(kx) - Sx*A(kz)
      Ay = A(ky) - Sy*A(kz)
      Bx = B(kx) - Sx*B(kz)
      By = B(ky) - Sy*B(kz)
      Cx = C(kx) - Sx*C(kz)
      Cy = C(ky) - Sy*C(kz)
      if(ldebug) print *,'local Triangle coords A:', Ax, Ay
      if(ldebug) print *,'local Triangle coords B:', Bx, By
      if(ldebug) print *,'local Triangle coords C:', Cx, Cy

      ! calculate scaled barycentric coordinates
      U = Cx*By - Cy*Bx;
      V = Ax*Cy - Ay*Cx;
      W = Bx*Ay - By*Ax;

      if(ldebug) print *,'Barycentric coords:', U, V, W

      ! fall back to test against edges using double precision
      if(ireal_dp.lt.ireal_dp) then
        if (any(approx([U,V,W],zero))) then
          CxBy = real(Cx, kind=ireal_dp) * real(By, kind=ireal_dp)
          CyBx = real(Cy, kind=ireal_dp) * real(Bx, kind=ireal_dp)
          U = real(CxBy - CyBx, kind=ireal_dp)
          AxCy = real(Ax, kind=ireal_dp) * real(Cy, kind=ireal_dp)
          AyCx = real(Ay, kind=ireal_dp) * real(Cx, kind=ireal_dp)
          V = real(AxCy - AyCx, kind=ireal_dp)
          BxAy = real(Bx, kind=ireal_dp) * real(Ay, kind=ireal_dp)
          ByAx = real(By, kind=ireal_dp) * real(Ax, kind=ireal_dp)
          W = real(BxAy - ByAx, kind=ireal_dp)
        endif
      endif

      !Perform edge tests. Moving this test before and at the end of the previous conditional gives higher performance.
      if(BACKFACE_CULLING) then
        if (U < zero .or. V < zero .or. W < zero) lhit=.False.
      else
        if ((U < zero .or. V < zero .or. W < zero) .and. &
          (U > zero .or. V > zero .or. W > zero)) lhit=.False.
      endif

      ! calculate determinant
      det = U + V + W
      if (.not.HIT_EDGE .and. approx(det, zero)) then
        if(ldebug) print *,'determinant zero: on edge?', det
        lhit=.False.
      endif

      !Calculate scaled z−coordinates of vertices and use them to calculate the hit distance.
      Az = Sz * A(kz)
      Bz = Sz * B(kz)
      Cz = Sz * C(kz)
      T = U * Az + V * Bz + W * Cz

      if(BACKFACE_CULLING) then
        if (T < zero .or. T > hit(4) * det) then
          if(ldebug) print *,'BACKFACE_CULLING T<0', T
          lhit = .False.
        endif
      else
        if(hit(4).ge.huge(hit)) then
          dist_times_det = sign(huge(det), det)
        else
          dist_times_det = hit(4)*det
        endif
        if(det < zero .and. ((T >= zero) .or. (T < dist_times_det))) then
          if(ldebug) print *,'det<0 && T>0', det, T
          lhit = .False.
        else if(det > zero .and. ((T <= zero) .or. (T > dist_times_det))) then
          if(ldebug) print *,'det>0 && T<0', det, T
          lhit = .False.
        endif
      endif

      ! normalize U, V, W, and T
      rcpDet = one / det
      b0 = U * rcpDet
      b1 = V * rcpDet
      b2 = W * rcpDet

      hit(1:3) = b0*tA + b1*tB + b2*tC
      hit(4) = T * rcpDet
      if(ldebug) print *,'Hit triangle', lhit, '::', hit
    end subroutine

      !> @brief determine if point is inside a rectangle p1,p2,p3
      function pnt_in_rectangle(p1,p2,p3, p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2,p3, p
        logical :: pnt_in_rectangle
        real(ireal_dp),parameter :: eps = epsilon(eps), eps2 = sqrt(eps)

        ! check for rectangular bounding box
        if ( p(1).lt.minval([p1(1),p2(1),p3(1)])-eps2 .or. p(1).gt.maxval([p1(1),p2(1),p3(1)])+eps2 ) then ! outside of xrange
            pnt_in_rectangle=.False.
            return
        endif
        if ( p(2).lt.minval([p1(2),p2(2),p3(2)])-eps2 .or. p(2).gt.maxval([p1(2),p2(2),p3(2)])+eps2 ) then ! outside of yrange
            pnt_in_rectangle=.False.
            return
        endif
        pnt_in_rectangle=.True.
      end function

      !> @brief determine if point is inside a triangle p1,p2,p3
      function pnt_in_triangle(p1,p2,p3, p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2,p3, p
        logical :: pnt_in_triangle
        real(ireal_dp),parameter :: eps = epsilon(eps)*10
        real(ireal_dp) :: a, b, c, edge_dist
        logical, parameter :: ldebug=.False.

        pnt_in_triangle = pnt_in_rectangle(p1,p2,p3, p)
        if(ldebug) print *,'pnt_in_triangle::pnt in rectangle:', p1, p2, p3, 'p', p, '::', pnt_in_triangle
        if (.not.pnt_in_triangle) then ! if pnt is not in rectangle, it is not in triangle!
          ! Then check for sides
          a = ((p2(2)- p3(2))*(p(1) - p3(1)) + (p3(1) - p2(1))*(p(2) - p3(2))) / ((p2(2) - p3(2))*(p1(1) - p3(1)) + (p3(1) - p2(1))*(p1(2) - p3(2)))
          b = ((p3(2) - p1(2))*(p(1) - p3(1)) + (p1(1) - p3(1))*(p(2) - p3(2))) / ((p2(2) - p3(2))*(p1(1) - p3(1)) + (p3(1) - p2(1))*(p1(2) - p3(2)))
          c = one - (a + b)

          pnt_in_triangle = all([a,b,c].ge.zero)
          if(ldebug) print *,'pnt_in_triangle::1st check:', a, b, c, '::', pnt_in_triangle
        endif

        if(.not.pnt_in_triangle) then
          pnt_in_triangle = pnt_in_triangle_convex_hull(p1,p2,p3, p)
          if(ldebug) print *,'pnt_in_triangle::convex hull:', pnt_in_triangle
        endif

        if(.not.pnt_in_triangle) then ! Compute distances to each edge and allow the check to be positive if the distance is small
          edge_dist = minval(distances_to_triangle_edges(p1,p2,p3,p))
          if(edge_dist.le.eps) then
            if((p(1).lt.min(p1(1),p2(1))) .or. (p(1).gt.max(p1(1),p2(1)) )) then
              ! is on line but ouside of segment
              continue
            else
              pnt_in_triangle=.True.
            endif
            if(ldebug) print *,'pnt_in_triangle edgedist:',edge_dist,'=>', pnt_in_triangle
          endif
        endif

        if(ldebug.and..not.pnt_in_triangle) print *,'pnt_in_triangle final:', pnt_in_triangle,'::',a,b,c,':',p, &
          'edgedist',distances_to_triangle_edges(p1,p2,p3,p),distances_to_triangle_edges(p1,p2,p3,p).le.eps
      end function

      function pnt_in_triangle_convex_hull(p1,p2,p3, p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2,p3, p
        logical :: pnt_in_triangle_convex_hull
        real(ireal_dp), dimension(2) :: v0, v1, v2
        real(ireal_dp) :: a,b

        v0 = p1
        v1 = p2-p1
        v2 = p3-p1

        a =  (cross_2d(p, v2) - cross_2d(v0, v2)) / cross_2d(v1, v2)
        b = -(cross_2d(p, v1) - cross_2d(v0, v1)) / cross_2d(v1, v2)

        pnt_in_triangle_convex_hull = all([a,b].ge.zero) .and. (a+b).le.one

        !print *,'points',p1,p2,p3,'::',p
        !print *,'a,b',a,b,'::',a+b, '::>',pnt_in_triangle_convex_hull
      end function

    pure function determine_normal_direction(normal, center_face, center_cell)
      ! return 1 if normal is pointing towards cell_center, -1 if its pointing
      ! away from it
      real(ireal_dp), intent(in) :: normal(:), center_face(:), center_cell(:)
      integer(iintegers) :: determine_normal_direction
      real(ireal_dp) :: dot
      dot = dot_product(normal, center_cell - center_face)
      determine_normal_direction = int(sign(one, dot), kind=iintegers)
    end function

      pure function distances_to_triangle_edges(p1,p2,p3,p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2,p3, p
        real(ireal_dp) :: distances_to_triangle_edges(3)
        distances_to_triangle_edges(1) = distance_to_edge(p1,p2,p)
        distances_to_triangle_edges(2) = distance_to_edge(p2,p3,p)
        distances_to_triangle_edges(3) = distance_to_edge(p1,p3,p)
      end function

      pure function distance_to_edge(p1,p2,p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2, p
        real(ireal_dp) :: distance_to_edge

        distance_to_edge = abs( (p2(2)-p1(2))*p(1) - (p2(1)-p1(1))*p(2) + p2(1)*p1(2) - p2(2)*p1(1) ) / norm(p2-p1)
      end function

      pure function rotate_angle_x(v,angle)
        ! left hand rule
        real(ireal_dp) :: rotate_angle_x(3)
        real(ireal_dp),intent(in) :: v(3), angle
        real(ireal_dp) :: M(3,3),s,c
        s=sin(deg2rad(angle))
        c=cos(deg2rad(angle))

        M(1,:)=[one ,zero ,zero]
        M(2,:)=[zero, c   , s  ]
        M(3,:)=[zero,-s   , c  ]

        rotate_angle_x = matmul(M,v)
      end function
      pure function rotate_angle_y(v,angle)
        ! left hand rule
        real(ireal_dp) :: rotate_angle_y(3)
        real(ireal_dp),intent(in) :: v(3), angle
        real(ireal_dp) :: M(3,3),s,c
        s=sin(deg2rad(angle))
        c=cos(deg2rad(angle))

        M(1,:)=[ c  ,zero , -s ]
        M(2,:)=[zero, one ,zero]
        M(3,:)=[ s  , zero, c  ]

        rotate_angle_y = matmul(M,v)
      end function
      pure function rotate_angle_z(v,angle)
        ! left hand rule
        real(ireal_dp) :: rotate_angle_z(3)
        real(ireal_dp),intent(in) :: v(3), angle
        real(ireal_dp) :: M(3,3),s,c
        s=sin(deg2rad(angle))
        c=cos(deg2rad(angle))

        M(1,:)=[ c  , s   ,zero]
        M(2,:)=[-s  , c   ,zero]
        M(3,:)=[zero, zero, one]

        rotate_angle_z = matmul(M,v)
      end function
    end module
