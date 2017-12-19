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

module f2c_tenstream

      use iso_c_binding

      use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,mpierr, zero, i0

      !use m_tenstream, only : init_tenstream, set_global_optical_properties, solve_tenstream, destroy_tenstream,&
      !                      tenstream_get_result, getvecpointer, restorevecpointer, &
      !                      t_coord,C_dir,C_diff,C_one,C_one_atm, C_one_atm1

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast,mean, CHKERR

      use m_pprts, only : init_pprts, t_solver, t_solver_8_10, t_solver_3_6, t_solver_1_2, &
        set_global_optical_properties, solve_pprts, destroy_pprts, &
        pprts_get_result_toZero, t_coord, petscVecToF90, petscGlobalVecToZero, f90VecToPetsc

#include "petsc/finclude/petsc.h"
      use petsc

      implicit none

      integer(mpiint) :: ierr
      class(t_solver), allocatable :: solver

contains

      subroutine tenstr_f2c_init(comm, solver_id, Nz,Nx,Ny,dx,dy,hhl, phi0, theta0, collapseindex) bind(C)
        ! initialize tenstream environment
        ! all nodes in communicator have to call this
        ! but only the zeroth node has to have meaningful values for the arguments except the communicator
        ! all but hhl is overwritten on nonzero nodes

        integer(c_int), value :: comm, solver_id
        integer(c_int),intent(inout) :: Nx,Ny,Nz
        real(c_double),intent(inout) :: dx,dy
        real(c_float), intent(in),dimension(Nz+1) :: hhl
        real(c_float), intent(inout) :: phi0,theta0
        integer(c_int),intent(inout) :: collapseindex

        integer(iintegers) :: oNx,oNy,oNz,ocollapseindex
        real(ireals) :: odx,ody,ophi0,otheta0
        real(ireals),allocatable :: ohhl(:)

        real(ireals),allocatable :: odz(:)
        integer(iintegers) :: k
        integer(mpiint) :: myid

        logical,save :: initialized=.False.

        if(initialized) return

        call init_mpi_data_parameters(comm)
        call read_commandline_options()
        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

        if(myid.eq.0) then
          oNx     = Nx
          oNy     = Ny
          oNz     = Nz
          odx     = dx
          ody     = dy
          ophi0   = phi0
          otheta0 = theta0
          ocollapseindex = collapseindex

          allocate( ohhl(size(hhl)) )
          ohhl = hhl
        endif

        call imp_bcast(comm, oNx    ,0_mpiint)
        call imp_bcast(comm, oNy    ,0_mpiint)
        call imp_bcast(comm, oNz    ,0_mpiint)
        call imp_bcast(comm, odx    ,0_mpiint)
        call imp_bcast(comm, ody    ,0_mpiint)
        call imp_bcast(comm, ophi0  ,0_mpiint)
        call imp_bcast(comm, otheta0,0_mpiint)
        call imp_bcast(comm, ocollapseindex,0_mpiint)

        call imp_bcast(comm, ohhl,0_mpiint)

        ! and overwrite input values to propagate values back to caller...
        Nx     = oNx
        Ny     = oNy
        Nz     = oNz
        dx     = odx
        dy     = ody
        phi0   = ophi0
        theta0 = otheta0
        collapseindex=ocollapseindex

        ! Now every process has the correct values
        print *,myid,'Initializing Tenstream environment from C Language :: domainshape',solver_id
        print *,myid,'Initializing Tenstream environment from C Language :: domainshape',oNx,oNy,oNz,'::',shape(ohhl)

        allocate(odz(oNz))
        do k=1,Nz
          odz(k) = ohhl(k) - ohhl(k+1)
        enddo

        select case(solver_id)
          case(0)
            allocate(t_solver_8_10::solver)
          case(1)
            allocate(t_solver_3_6::solver)
          case(2)
            allocate(t_solver_1_2::solver)
        end select

        call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, solver, dz1d=odz)

        initialized=.True.
      end subroutine

      subroutine tenstr_f2c_set_global_optical_properties(Nz, Nx, Ny, albedo, kabs, ksca, g, planck) bind(c)
        integer(c_int), value :: Nx,Ny,Nz
        real(c_float), value :: albedo
        real(c_float),intent(in),dimension(Nz  ,Nx,Ny) :: kabs, ksca, g
        real(c_float),intent(in),dimension(Nz+1,Nx,Ny) :: planck

        real(ireals) :: oalbedo
        real(ireals),allocatable,dimension(:,:,:) :: okabs, oksca, og, oplanck

        oalbedo = albedo

        if(solver%myid.eq.0) then
          allocate( okabs  (Nz  ,Nx,Ny) ); okabs   = kabs
          allocate( oksca  (Nz  ,Nx,Ny) ); oksca   = ksca
          allocate( og     (Nz  ,Nx,Ny) ); og      = g
          allocate( oplanck(Nz+1,Nx,Ny) ); oplanck = planck

          if(any(oplanck.gt.zero)) then
            call set_global_optical_properties(solver, oalbedo, okabs, oksca, og, oplanck)
          else
            call set_global_optical_properties(solver, oalbedo, okabs, oksca, og)
          endif

          print *,'mean kabs  ',sum(okabs)  /size(okabs)
          print *,'mean ksca  ',sum(oksca)  /size(oksca)
          print *,'mean g     ',sum(og)     /size(og)
          print *,'mean planck',sum(oplanck)/size(oplanck)
        else !slave
          call set_global_optical_properties(solver, oalbedo)
        endif
      end subroutine

      subroutine tenstr_f2c_solve(comm, edirTOA) bind(c)
        ! solve tenstream equations
        ! optical properties have had to be set and environment had to be initialized
        ! incoming solar radiation need only be set by zeroth node
        integer(c_int), value :: comm
        real(c_float), value :: edirTOA
        real(ireals) :: oedirTOA

        if(solver%myid.eq.0) oedirTOA = edirTOA
        call imp_bcast(comm, oedirTOA, 0_mpiint)

        call solve_pprts(solver, oedirTOA)
      end subroutine

      subroutine tenstr_f2c_destroy() bind(c)
        call destroy_pprts(solver, lfinalizepetsc=.False.)
        deallocate(solver)
      end subroutine

      subroutine tenstr_f2c_get_result(Nz,Nx,Ny, res_edir,res_edn,res_eup,res_abso) bind(c)
        ! after solving equations -- retrieve the results for edir,edn,eup and absorption
        ! only zeroth node gets the results back.

        integer(c_int), value :: Nx,Ny,Nz
        real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edir
        real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edn
        real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_eup
        real(c_float),intent(out),dimension(Nz  ,Nx,Ny) :: res_abso
        real(ireals),allocatable,dimension(:,:,:) :: res

        real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso

        call pprts_get_result_toZero(solver,redir,redn,reup,rabso)

     end subroutine
end module
