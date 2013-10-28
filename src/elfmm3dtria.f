cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains the FMM routines for elastostatic layer
c        potentials in free space in R^3. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      elfmm3dtria - evaluates the elastostatic potential ON SURFACE due
c         to a collection of flat triangles with piecewise constant
c         single and/or double layer densities using the Fast Multipole
c         Method.
c
c      elfmm3dtriatarg - evaluates the elastostatic potential ON OR OFF
c         SURFACE due to a collection of flat triangles with constant
c         single and/or double layer densities using the Fast Multipole
c         Method.
c
c      el3dtriadirecttarg - evaluates the elastostatic potential ON OR
c         OFF SURFACE due to a collection of flat triangles with
c         constant single and/or double layer densities using the direct
c         O(N^2) algorithm.
c
c
c      elust3triadirectself - evaluates the elastostatic potential at a
c         single surface point (triangle centroid) due to a collection
c         of flat triangles with piecewise constant single layer density
c         BY DIRECT CALCULATION.
c
c      elust3triadirecttarg - evaluates the elastostatic potential at an
c         (OFF SURFACE) target due to a collection of flat triangles
c         with piecewise constant single layer density BY DIRECT
c         CALCULATION.
c
c      eltst3triadirectself - evaluates the elastostatic potential at a
c         single surface point (triangle centroid) due to a collection
c         of flat triangles with piecewise constant double layer density
c         BY DIRECT CALCULATION.
c     
c      eltst3triadirecttarg - evaluates the elastostatic potential at an
c         (OFF SURFACE) target due to a collection of flat triangles
c         with piecewise constant double layer density BY DIRECT
c         CALCULATION.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine elfmm3dtria
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain)
c
c     This subroutine evaluates the elastostatic potential 
c     due to a collection of flat triangles with constant
c     single and/or double layer densities using the Fast Multipole Method.
c
c     Unfortunately, in different communities the source strengths
c     are referred to by different names. Here we will use the 
c     mathematical conventions of single and double layer potentials.
c     
c   
c     The single layer Green's function maps 
c     traction (surface force density) to displacement/strain.
c     (In some classical literature, this is called the <<single force>>.)
c     We will use SIGMA_SL to describe the source density.
c
c     (The resulting displacement is continuous across the surface.)
c
c     The double layer Green's function maps 
c     generalized slip/jump in displacement to displacement/strain
c     (In some classical literature, this is called the <<double force>>.)
c     We will use SIGMA_DL to describe the source density.
c
c     (The resulting displacement is NOT continuous across the surface.)
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of single layer source 
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of double layer source
c
c     iprec:  FMM precision flag
c
c             -2 => tolerance =.5d0
c             -1 => tolerance =.5d-1
c             0 => tolerance =.5d-2
c             1 => tolerance =.5d-3
c             2 => tolerance =.5d-6
c             3 => tolerance =.5d-9
c             4 => tolerance =.5d-12
c             5 => tolerance =.5d-15
c
c     These errors were established for the electrostatic case and are
c     approximately valid for displacement. For strain, there is
c     approximately one digit of loss in the above table.
c
c     OUTPUT:
c
c     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
c
c     The main FMM routine permits both evaluation on surface
c     and at a collection of off-surface targets. 
c     This subroutine is used to simplify the user interface 
c     (by setting the number of targets to zero) and calling the more 
c     general FMM.
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        integer nparts,ntargs

        ntargs=0
        ifptfrctarg=0
        ifstraintarg=0
        call elfmm3dtriatarg
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)

        return
        end
c
c
c
c
c
      subroutine elfmm3dtriatarg
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
c
c       This is the principal subroutine for evaluating elastostatic
c       layer potentials on (flat) triangulated surfaces.  It permits
c       the evaluation of a single layer potential with piecewise
c       constant density defined by the (force) vector sigma_sl and a
c       dipole layer potential with piecewise constant density vector
c       sigma_dl and dipole orientation defined by trinorm.
c
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutines 
c       elfmm3dtriatargmain_fast  for the far field interactions and 
c       and elfmm3dtriatarg0  for the near field interactions.
c 
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
c     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrctarg(3,ntargs),straintarg(3,3,ntargs)
        integer nparts,ntargs
        real *8, allocatable :: w(:)
c       
        ier=0
        lused=0
c
c       ... allocate work arrays
c
        icharge=lused+1
        lcharge=2*nparts *3
        lused=lused+lcharge

        idipstr=lused+1
        ldipstr=2*nparts *3
        lused=lused+ldipstr

        idipvec=lused+1
        ldipvec=3*nparts *3
        lused=lused+ldipvec

        ipot=lused+1
        lpot=2*nparts
        lused=lused+lpot

        ifld=lused+1
        lfld=2*3*nparts
        lused=lused+lfld

        ihess=lused+1
        lhess=2*6*nparts
        lused=lused+lhess

        ihessmatr=lused+1
        lhessmatr=3*3*nparts
        lused=lused+lhessmatr

        ipottarg=lused+1
        lpottarg=2*ntargs
        lused=lused+lpottarg

        ifldtarg=lused+1
        lfldtarg=2*3*ntargs
        lused=lused+lfldtarg

        ihesstarg=lused+1
        lhesstarg=2*6*ntargs
        lused=lused+lhesstarg

        ihessmatrtarg=lused+1
        lhessmatrtarg=3*3*ntargs
        lused=lused+lhessmatrtarg
c
        allocate( w(lused), stat=ier)
        if( ier .ne. 0 ) return
c
c     call FMM to account for all far field interactions.
c     
c     The subroutine elfmm3dtriatargmain_fast makes 4 calls
c     to a scalar (Laplace FMM) according to the canonical
c     decomposition of the elastostatic Green fucntions.
c
c
        call elfmm3dtriatargmain_fast
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,trinorm,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg,
     $     w(icharge),w(idipstr),w(idipvec),
     $     w(ipot),w(ifld),w(ihess),w(ihessmatr),
     $     w(ipottarg),w(ifldtarg),w(ihesstarg),w(ihessmatrtarg))
c
c     reconstruct FMM data structure and account for  all local 
c     interactions using quadrature routines for piecewise
c     constant densities - no interactions are saved in the 
c     present version.
c     
        call elfmm3dtriatarg0
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,trinorm,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
c
        return
        end
c
c
c
c
c
C*********************************
      subroutine elfmm3dtriatargmain_fast
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg,
     $     charge,dipstr,dipvec,pot,fld,hess,hessmatr,
     $     pottarg,fldtarg,hesstarg,hessmatrtarg)
c
c     FMM calculation subroutine for elastostatic N-body problem
c
c     4 Laplace FMM calls, uses linear charge and dipole densities
c
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c             -2 => tolerance =.5d0
c             -1 => tolerance =.5d-1
c             0 => tolerance =.5d-2
c             1 => tolerance =.5d-3
c             2 => tolerance =.5d-6
c             3 => tolerance =.5d-9
c             4 => tolerance =.5d-12
c             5 => tolerance =.5d-15
c
c
c     OUTPUT:
c
c     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
c     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,nparts),trinorm(3,nparts)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrctarg(3,ntargs),straintarg(3,3,ntargs)
        integer nparts,ntargs
c       
        complex *16 charge(3,1)
        complex *16 dipstr(3,1)
        real *8 dipvec(3,3,1)
c
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1) 
        real *8 hessmatr(3,3,1)

        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        complex *16 hesstarg(6,1) 
        real *8 hessmatrtarg(3,3,1)
c
c
        do k=1,nparts
c
        if( ifptfrc .eq. 1 ) then
           ptfrc(1,k) = 0.0d0
           ptfrc(2,k) = 0.0d0
           ptfrc(3,k) = 0.0d0
        endif
C
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
           hessmatr(i,j,k) = 0.0d0
        enddo
        enddo
        endif        
c
        enddo
c
c
        do k=1,ntargs
c
        if( ifptfrctarg .eq. 1 ) then
           ptfrctarg(1,k) = 0.0d0
           ptfrctarg(2,k) = 0.0d0
           ptfrctarg(3,k) = 0.0d0
        endif
C
        if( ifstraintarg .eq. 1 ) then
        do i=1,3
        do j=1,3
           hessmatrtarg(i,j,k) = 0.0d0
        enddo
        enddo
        endif
c
        enddo
c
c
c
        ifpot=0
        iffld=0
        ifhess=0
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
        ifpot=1
        iffld=1
        endif
        if( ifstrain .eq. 1 ) ifhess=1

        ifpottarg=0
        iffldtarg=0
        ifhesstarg=0
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
        ifpottarg=1
        iffldtarg=1
        endif
        if( ifstraintarg .eq. 1 ) ifhesstarg=1
c
c
        c1 = (rlam+rmu)/(rlam+2*rmu)
        c2 = (rmu)/(rlam+2*rmu)
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
c       Note, that the constructed charge and dipole vectors are constant
c       functions on the triangle,  we call lfmm3dtrilhesstarg routine 
c       to simplify the things (we can use lfmm3dtrilhesstarg here)
c       
        do j = 1,3

           ifcharge=0
           ifdipole=0

           do m=1,3
           do k = 1,nparts
              charge(m,k) = 0
              dipstr(m,k) = 0
              dipvec(1,m,k) = 0
              dipvec(2,m,k) = 0
              dipvec(3,m,k) = 0
              if( ifsingle .eq. 1 ) then
                charge(m,k) = sigma_sl(j,k)/(2*rmu)
                ifcharge=1
              endif
              if( ifdouble .eq. 1 ) then
                DIPSTR(m,k) = 1
                dipvec(1,m,k)=sigma_dv(1,k)*sigma_dl(j,k)
                dipvec(2,m,k)=sigma_dv(2,k)*sigma_dl(j,k)
                dipvec(3,m,k)=sigma_dv(3,k)*sigma_dl(j,k)
                dipvec(1,m,k)=dipvec(1,m,k)+sigma_dl(1,k)*sigma_dv(j,k)
                dipvec(2,m,k)=dipvec(2,m,k)+sigma_dl(2,k)*sigma_dv(j,k)
                dipvec(3,m,k)=dipvec(3,m,k)+sigma_dl(3,k)*sigma_dv(j,k)
                dipvec(1,m,k)=dipvec(1,m,k)/2
                dipvec(2,m,k)=dipvec(2,m,k)/2
                dipvec(3,m,k)=dipvec(3,m,k)/2
                ifdipole=1
              endif
           enddo
           enddo
c
           call lfmm3dtrilhesstarg(ier,iprec,
     $        nparts,triangle,trinorm,source,
     $        ifcharge,charge,ifdipole,dipstr,dipvec,
     $        ifpot,pot,iffld,fld,ifhess,hess,
     $        ntargs,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $        ifhesstarg,hesstarg)
c
           call elfmm3dlap1(nparts,j,c1,c2,pot,fld,hess,
     $        source,ifptfrc,ptfrc,ifstrain,hessmatr)
           call elfmm3dlap1(ntargs,j,c1,c2,pottarg,fldtarg,hesstarg,
     $        target,ifptfrctarg,ptfrctarg,ifstraintarg,hessmatrtarg)
        enddo
c
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c
c       Note, that the constructed charge and dipole vectors are linear
c       functions on the triangle, we must call lfmm3dtrilhesstarg routine
c       here (a significant loss of accuracy in near field will occur if
c       we use lfmm3dtrilhesstarg instead)
c       
        ifcharge=0
        ifdipole=0
c        
        do m=1,3
        do k = 1,nparts
          charge(m,k) = 0
          dipstr(m,k) = 0
          dipvec(1,m,k) = 0
          dipvec(2,m,k) = 0
          dipvec(3,m,k) = 0
          if( ifsingle .eq. 1 ) then
          charge(m,k) = 
     $      (sigma_sl(1,k)*triangle(1,m,k)+
     $       sigma_sl(2,k)*triangle(2,m,k)+
     $       sigma_sl(3,k)*triangle(3,m,k))*c1/(2*rmu)
          ifcharge = 1
          endif
          if( ifdouble .eq. 1 ) then
          charge(m,k) = charge(m,k) + 
     $        (sigma_dl(1,k)*sigma_dv(1,k) + 
     1         sigma_dl(2,k)*sigma_dv(2,k) + 
     2         sigma_dl(3,k)*sigma_dv(3,k))*c2
          dipstr(m,k) = 1 
          dipvec(1,m,k) = c1*sigma_dv(1,k)*
     $        (sigma_dl(1,k)*triangle(1,m,k) + 
     1         sigma_dl(2,k)*triangle(2,m,k) + 
     2         sigma_dl(3,k)*triangle(3,m,k) )
          dipvec(2,m,k) = c1*sigma_dv(2,k)*
     $        (sigma_dl(1,k)*triangle(1,m,k) + 
     1         sigma_dl(2,k)*triangle(2,m,k) + 
     2         sigma_dl(3,k)*triangle(3,m,k) )
          dipvec(3,m,k) = c1*sigma_dv(3,k)*
     $        (sigma_dl(1,k)*triangle(1,m,k) + 
     1         sigma_dl(2,k)*triangle(2,m,k) + 
     2         sigma_dl(3,k)*triangle(3,m,k) )
          dipvec(1,m,k) = dipvec(1,m,k) + c1*sigma_dl(1,k)*
     $        (sigma_dv(1,k)*triangle(1,m,k) + 
     1         sigma_dv(2,k)*triangle(2,m,k) + 
     2         sigma_dv(3,k)*triangle(3,m,k))
          dipvec(2,m,k) = dipvec(2,m,k) + c1*sigma_dl(2,k)*
     $        (sigma_dv(1,k)*triangle(1,m,k) + 
     1         sigma_dv(2,k)*triangle(2,m,k) + 
     2         sigma_dv(3,k)*triangle(3,m,k))
          dipvec(3,m,k) = dipvec(3,m,k) + c1*sigma_dl(3,k)*
     $        (sigma_dv(1,k)*triangle(1,m,k) + 
     1         sigma_dv(2,k)*triangle(2,m,k) + 
     2         sigma_dv(3,k)*triangle(3,m,k))
          dipvec(1,m,k) = dipvec(1,m,k)/2
          dipvec(2,m,k) = dipvec(2,m,k)/2
          dipvec(3,m,k) = dipvec(3,m,k)/2
          ifcharge = 1
          ifdipole = 1          
          endif
        enddo
        enddo
c
        call lfmm3dtrilhesstarg(ier,iprec,
     $     nparts,triangle,trinorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntargs,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c
        call elfmm3dlap2(nparts,pot,fld,hess,
     $     ifptfrc,ptfrc,ifstrain,hessmatr)
        call elfmm3dlap2(ntargs,pottarg,fldtarg,hesstarg,
     $     ifptfrctarg,ptfrctarg,ifstraintarg,hessmatrtarg)
c
c
 2000   continue
c
        do k=1,nparts
           if( ifptfrc .eq. 1 ) then
              ptfrc(1,k)=ptfrc(1,k)
              ptfrc(2,k)=ptfrc(2,k)
              ptfrc(3,k)=ptfrc(3,k)
           endif
c       
           if( ifstrain .eq. 1 ) then
              do i=1,3
              do j=1,3
                 strain(i,j,k)=(hessmatr(i,j,k)+hessmatr(j,i,k))/2
              enddo
              enddo
           endif
        enddo
c 
        do k=1,ntargs
c       
        if( ifptfrctarg .eq. 1 ) then
           ptfrctarg(1,k)=ptfrctarg(1,k)
           ptfrctarg(2,k)=ptfrctarg(2,k)
           ptfrctarg(3,k)=ptfrctarg(3,k)
        endif
c       
        if( ifstraintarg .eq. 1 ) then
           do i=1,3
           do j=1,3
              straintarg(i,j,k)=
     $           (hessmatrtarg(i,j,k)+hessmatrtarg(j,i,k))/2
           enddo
           enddo
        endif
c
        enddo
c       
        return
        end
C
c
c
c
c
c
c
c
        subroutine elfmm3dtriatarg0(ier,iprec,rlam,rmu,
     $     triangle,trinorm,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,
     $     ifptfrctarg,ptfrctarg,ifstraintarg,straintarg)
c
c     FMM calculation subroutine for elastostatic N-body problem
c
c     Direct evaluation routine, for local interactions only
c
c     Constant densities on flat triangles.
c     Note that currently, SIGMA_DV must be the same as TRINORM.
c
C     INPUT:
C
C     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     nsource = number of sources
C     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
C     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
C     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
C     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
C     target(3,ntargs) = evaluation target points
C
c     iprec:  FMM precision flag
c
C     OUTPUT:
C
C     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
C     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
c
        implicit real *8 (a-h,o-z)
        dimension triangle(3,3,1),trinorm(3,1)
        dimension source(3,1)
        dimension sigma_sl(3,1)
        dimension sigma_dl(3,1)
        dimension sigma_dv(3,1)
        dimension ptfrc(3,1)
        dimension strain(3,3,1)
        dimension target(3,1)
        dimension ptfrctarg(3,1)
        dimension straintarg(3,3,1)
c
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
c
        dimension timeinfo(10)
c       
        dimension center(3)
c       
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
c       
        integer box(20)
        dimension center0(3),corners0(3,8)
c       
        integer box1(20)
        dimension center1(3),corners1(3,8)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
c
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c       ... build the oct-tree
c       
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c       
        call prin2('epsfmm=*',epsfmm,1)
c
        if( iprec .eq. -2 ) nbox=8/3
        if( iprec .eq. -1 ) nbox=15/3
        if( iprec .eq. 0 ) nbox=30/3
        if( iprec .eq. 1 ) nbox=60/3
        if( iprec .eq. 2 ) nbox=120/3
        if( iprec .eq. 3 ) nbox=240/3
        if( iprec .eq. 4 ) nbox=480/3
        if( iprec .eq. 5 ) nbox=700/3
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        call prinf('nbox=*',nbox,1)
c
c
c     create oct-tree data structure
c
        ntot = 100*(nsource+ntarget)+10000
        if( iprec .eq. -2 ) ntot = ntot * 1.5*1.5*1.5
        if( iprec .eq. -1 ) ntot = ntot * 1.5*1.5
        do ii = 1,10
           allocate (wlists(ntot))
           call lfmm3dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
              goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return          
        endif
c
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
c       
        call prin2('scale=*',scale,nlev+1)
c       
c       
c       carve up workspace further
c
c     itriaflatsort is pointer for sorted triangle coordinates
c     itrianormsort is pointer for sorted triangle normals
c     isourcesort is pointer for sorted triangle centroids
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
c
        itrianglesort = lused7 
        ltrianglesort = 3*3*nsource
        itrinormsort = itrianglesort + ltrianglesort
        ltrinormsort = 3*nsource

        isourcesort = itrinormsort + ltrinormsort 
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget

        isigma_slsort = itargetsort+ltargetsort
        if (ifsingle.eq.1) then
          lsigma_slsort = 3*nsource
        else
          lsigma_slsort = 3
        endif
        isigma_dlsort = isigma_slsort+lsigma_slsort
        if (ifdouble.eq.1) then
          lsigma_dlsort = 3*nsource
        else
          lsigma_dlsort = 3
        endif
        isigma_dvsort = isigma_dlsort+lsigma_dlsort
        if (ifdouble.eq.1) then
          lsigma_dvsort = 3*nsource
        else
          lsigma_dvsort = 3
        endif

        lused7 = isigma_dvsort+lsigma_dvsort
c
c
c       ... allocate the potential and field arrays
c
c
        iptfrc = lused7 
        if( ifptfrc .eq. 1) then
        lptfrc = 2*(3*nsource)
        else
        lptfrc=6
        endif
        lused7=lused7+lptfrc
c      
        istrain = lused7
        if( ifstrain .eq. 1) then
        lstrain = 2*(3*3*nsource)
        else
        lstrain= 2*3*3
        endif
        lused7=lused7+lstrain
c      
        iptfrctarg = lused7
        if( ifptfrctarg .eq. 1) then
        lptfrctarg = 2*(3*ntarget)
        else
        lptfrctarg=6
        endif
        lused7=lused7+lptfrctarg
c      
        istraintarg = lused7
        if( ifstraintarg .eq. 1) then
        lstraintarg = 2*(3*3*ntarget)
        else
        lstraintarg= 2*3*3
        endif
        lused7=lused7+lstraintarg
c      
c
        call prinf(' lused7 is *',lused7,1)
c
c   
c       ... allocate temporary arrays
c
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                   lused7 is *',lused7,1)
           ier = 8
           return          
        endif
c

        call l3dreordertarg
     $     (nsource,source,wlists(iisource),w(isourcesort))
        if( ifsingle .eq. 1 ) then
        call l3dreordertarg
     $     (nsource,sigma_sl,wlists(iisource),w(isigma_slsort))
        endif
        if( ifdouble .eq. 1 ) then
        call l3dreordertarg
     $     (nsource,sigma_dl,wlists(iisource),w(isigma_dlsort))
        call l3dreordertarg
     $     (nsource,sigma_dv,wlists(iisource),w(isigma_dvsort))
        endif
        
        call l3dreordertria(nsource,wlists(iisource),
     $     triangle,w(itrianglesort),trinorm,w(itrinormsort))
c
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
c
c
c
        call elfmm3dtriatarg0_evalloc(ier,iprec,RLAM,RMU,
     $     w(itrianglesort),w(itrinormsort),
     $     nsource,w(isourcesort),
     $     ifsingle,w(isigma_slsort),
     $     ifdouble,w(isigma_dlsort),w(isigma_dvsort),
     $     ifptfrc,w(iptfrc),ifstrain,w(istrain),
     $     ntarget,w(itargetsort),
     $     ifptfrctarg,w(iptfrctarg),ifstraintarg,w(istraintarg),
     $     nboxes,laddr,nlev,wlists(iwlists),lwlists)
c
        call prinf('lwlists=*',lwlists,1)
        call prinf('lused total =*',lused7,1)
c       
        call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifptfrc .eq. 1) 
     $     call l3dptsort(nsource,wlists(iisource),w(iptfrc),ptfrc)
        if(ifstrain .eq. 1) 
     $     call l3dstsort(nsource,wlists(iisource),w(istrain),strain)
c
        if(ifptfrctarg .eq. 1 )
     $     call l3dptsort(ntarget,wlists(iitarget),
     $     w(iptfrctarg),ptfrctarg)
        if(ifstraintarg .eq. 1) 
     $     call l3dstsort(ntarget,wlists(iitarget),
     $     w(istraintarg),straintarg)
c       
        return
        end
c
c
c
c
c
        subroutine l3dptsort(n,isource,ptfrcsort,ptfrc)
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        real *8 ptfrc(3,1),ptfrcsort(3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        do m=1,3
ccc        ptfrc(m,isource(i))=ptfrcsort(m,i)
        ptfrc(m,isource(i))=ptfrc(m,isource(i))+ptfrcsort(m,i)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l3dstsort(n,isource,strainsort,strain)
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        real *8 strain(3,3,1),strainsort(3,3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        do j=1,3
        do m=1,3
ccc        strain(m,j,isource(i))=strainsort(m,j,i)
        strain(m,j,isource(i))=
     $     strain(m,j,isource(i))+strainsort(m,j,i)
        enddo
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine elfmm3dtriatarg0_evalloc(ier,iprec,rlam,rmu,
     $     triangle,trinorm,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,
     $     ifptfrctarg,ptfrctarg,ifstraintarg,straintarg,
     $     nboxes,laddr,nlev,wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension triangle(3,3,1),trinorm(3,1)
        dimension source(3,1)
        dimension sigma_sl(3,1)
        dimension sigma_dl(3,1)
        dimension sigma_dv(3,1)
        dimension ptfrc(3,1)
        dimension strain(3,3,1)
        dimension target(3,1)
        dimension ptfrctarg(3,1)
        dimension straintarg(3,3,1)
c
        dimension laddr(2,200)
        dimension list(10 000)
c
        dimension timeinfo(10)
c
        dimension wlists(1)
c
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
c
ccc        save
c
c     
c       ... set the displacement and strain to zero
c
        do i=1,nsource
        if (ifptfrc .eq. 1) then
        ptfrc(1,i)=0
        ptfrc(2,i)=0
        ptfrc(3,i)=0
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,i)=0
        strain(2,1,i)=0
        strain(3,1,i)=0
        strain(1,2,i)=0
        strain(2,2,i)=0
        strain(3,2,i)=0
        strain(1,3,i)=0
        strain(2,3,i)=0
        strain(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,i)=0
        ptfrctarg(2,i)=0
        ptfrctarg(3,i)=0
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,i)=0
        straintarg(2,1,i)=0
        straintarg(3,1,i)=0
        straintarg(1,2,i)=0
        straintarg(2,2,i)=0
        straintarg(3,2,i)=0
        straintarg(1,3,i)=0
        straintarg(2,3,i)=0
        straintarg(3,3,i)=0
        endif
        enddo
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
        call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        ifprint=0
        if (ifprint .eq. 1) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .eq. 1) then
               call prinf('npts=*',npts,1)
            endif
        endif
c
c
        if (nkids .eq. 0 ) then
c
c       ... evaluate self interactions
c
        call elfmm3dtria_direct_self(box,
     $     rlam,rmu,triangle,trinorm,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .eq. 1) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c    
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c
               call elfmm3dtria_direct(box1,box,
     $            rlam,rmu,triangle,trinorm,source,
     $            ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $            ifptfrc,ptfrc,ifstrain,strain,
     $            target,ifptfrctarg,ptfrctarg,
     $            ifstraintarg,straintarg)
c
 6203           continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
c
ccc        call prin2('inside fmm, ptfrc=*',ptfrc,3*nsource)
ccc        call prin2('inside fmm, ptfrctarg=*',ptfrctarg,3*ntarget)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        call prin2('timeinfo=*',timeinfo,8)
c       
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
c       
        return
        end
c
c
c
c
c
        subroutine elfmm3dtria_direct(box,box1,
     $     rlam,rmu,triangle,trinorm,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triangle(3,3,1),trinorm(3,1),source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
c
c       ... sources
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptfrc0,strain0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(14),box1(14)+box1(15)-1
        if (ifsingle .eq. 1 ) then
        call elust3triadirecttarg
     $     (rlam,rmu,box(15),triangle(1,1,box(14)),
     $     sigma_sl(1,box(14)),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        call eltst3triadirecttarg
     $     (rlam,rmu,box(15),triangle(1,1,box(14)),
     $     sigma_dl(1,box(14)),sigma_dv(1,box(14)),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptfrc0,strain0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(16),box1(16)+box1(17)-1
        if (ifsingle .eq. 1 ) then
        call elust3triadirecttarg
     $     (rlam,rmu,box(15),triangle(1,1,box(14)),
     $     sigma_sl(1,box(14)),
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        call eltst3triadirecttarg
     $     (rlam,rmu,box(15),triangle(1,1,box(14)),
     $     sigma_dl(1,box(14)),sigma_dv(1,box(14)),
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine elfmm3dtria_direct_self(box,
     $     rlam,rmu,triangle,trinorm,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triangle(3,3,1),trinorm(3,1),source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
c
c       ... sources
c
        ione=1
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptfrc0,strain0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(14),box(14)+box(15)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifsingle .eq. 1 ) then
        if( i .eq. j ) then
        call elust3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        else
        call elust3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        if( i .eq. j ) then 
        call eltst3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        else
        call eltst3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptfrc0,strain0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(16),box(16)+box(17)-1
        if (ifsingle .eq. 1 ) then
        call elust3triadirecttarg
     $     (rlam,rmu,box(15),triangle(1,1,box(14)),
     $     sigma_sl(1,box(14)),
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        call eltst3triadirecttarg
     $     (rlam,rmu,box(15),triangle(1,1,box(14)),
     $     sigma_dl(1,box(14)),sigma_dv(1,box(14)),
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine el3dtriadirecttarg(
     $     rlam,rmu,triangle,trinorm,nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain,ntarget,
     $     target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
        implicit real *8 (a-h,o-z)
c
c
c       Elastostatic interactions in R^3: evaluate all pairwise triangle
c       interactions and interactions with targets using the direct
c       O(N^2) algorithm.
c
c       INPUT:
c
c       rlam,rmu - Lame parameters
c       triangle(3,3,nsource) - array of triangles in standard format
c       trinorm(3,nsource) - array of triangle normals
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifptfrc - displacement computation flag
c       ifstrain - strain computation flag
c       ifptfrctarg - target displacement computation flag
c       ifstraintarg - target strain computation flag
c
c       OUTPUT:
c
c       ptfrc(3,nsource) - displacement at source locations
c       strain(3,3,nsource) - strain at source locations
c       ptfrctarg(3,ntarget) - displacement at target locations
c       straintarg(3,3,ntarget) - strain at target locations
c
c
        dimension triangle(3,3,1),trinorm(3,1),source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
c
c
c     NOTE: In the present version, dipole vectors SIGMA_DV must be SET EQUAL
c     to the triangle normal for elt3d direct routines
c
        do i=1,nsource
        if( ifptfrc .eq. 1) then
           ptfrc(1,i)=0
           ptfrc(2,i)=0
           ptfrc(3,i)=0
        endif
        if( ifstrain .eq. 1) then
           strain(1,1,i)=0
           strain(2,1,i)=0
           strain(3,1,i)=0
           strain(1,2,i)=0
           strain(2,2,i)=0
           strain(3,2,i)=0
           strain(1,3,i)=0
           strain(2,3,i)=0
           strain(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifptfrctarg .eq. 1) then
           ptfrctarg(1,i)=0
           ptfrctarg(2,i)=0
           ptfrctarg(3,i)=0
        endif
        if( ifstraintarg .eq. 1) then
           straintarg(1,1,i)=0
           straintarg(2,1,i)=0
           straintarg(3,1,i)=0
           straintarg(1,2,i)=0
           straintarg(2,2,i)=0
           straintarg(3,2,i)=0
           straintarg(1,3,i)=0
           straintarg(2,3,i)=0
           straintarg(3,3,i)=0
        endif
        enddo
c
c       ... sources
c
        ione=1
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do j=1,nsource
        do i=1,nsource
c
        if (ifsingle .eq. 1 ) then
        if( i .eq. j ) then
        call elust3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        else
        call elust3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        if( i .eq. j ) then 
        call eltst3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        else
        call eltst3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
        enddo
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do j=1,ntarget
        if (ifsingle .eq. 1 ) then
        call elust3triadirecttarg
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        call eltst3triadirecttarg
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
C$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
c***********************************************************************
c
c       Quadrature routines for elastostatic single and double layers
c       Constant densities on flat triangles
c
c***********************************************************************
c
c
        subroutine elust3triadirecttarg_one
     $     (rlam,rmu,triangle,sigma_sl,ifself,target,
     $     ptfrc,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c
c     Direct evaluation of displacement and strain due to constant
c     single layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET
c     due to piecewise-constant single layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3)          SLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        dimension triangle(3,3),sigma_sl(3)
        dimension ptfrc(3),strain(3,3)
        dimension w(20),vert1(3),vert2(3),vert3(3)
        dimension vectout(3),vertout(3), ders(3,3)
        dimension rtable(0:2,0:2,0:2),btable(0:4,0:4,0:4)
c
        c1 = (rlam+rmu)/(rlam+2*rmu)
        c2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
        call tri_for_vect(w,sigma_sl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0
c
        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif


        if( ifstrain .eq. 0 ) maxb=2
        if( ifstrain .eq. 0 ) maxr=0
        if( ifstrain .eq. 1 ) maxb=3
        if( ifstrain .eq. 1 ) maxr=1

        do k=0,maxb
        do i=0,maxb
        do j=0,maxb
        if( i+j+k .le. maxb ) then
        call triabtable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        btable(i,j,k)=d
        endif
        enddo
        enddo
        enddo
        
        do k=0,maxr
        do i=0,maxr
        do j=0,maxr
        if( i+j+k .le. maxr ) then
        call triartable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        rtable(i,j,k)=d
        endif
        enddo
        enddo
        enddo


        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(2,0,0)
        valx = -d*c1

        d=btable(1,1,0)
        valy = -d*c1

        d=btable(1,0,1)
        valz = -d*c1

        d=rtable(0,0,0)
        valx = valx+d*2

        if( ifstrain .eq. 1 ) then
        d=btable(3,0,0)
        valxx = -d*c1

        d=btable(2,1,0)
        valxy = -d*c1

        d=btable(2,0,1)
        valxz = -d*c1

        d=btable(2,1,0)
        valyx = -d*c1

        d=btable(1,2,0)
        valyy = -d*c1

        d=btable(1,1,1)
        valyz = -d*c1

        d=btable(2,0,1)
        valzx = -d*c1

        d=btable(1,1,1)
        valzy = -d*c1

        d=btable(1,0,2)
        valzz = -d*c1

        d=rtable(1,0,0)
        valxx = valxx+d*2

        d=rtable(0,1,0)
        valxy = valxy+d*2

        d=rtable(0,0,1)
        valxz = valxz+d*2
        endif

        endif

        if ( k .eq. 2 ) then
        
        d=btable(1,1,0)
        valx = -d*c1

        d=btable(0,2,0)
        valy = -d*c1

        d=btable(0,1,1)
        valz = -d*c1

        d=rtable(0,0,0)
        valy = valy+d*2

        if( ifstrain .eq. 1 ) then
        d=btable(2,1,0)
        valxx = -d*c1

        d=btable(1,2,0)
        valxy = -d*c1

        d=btable(1,1,1)
        valxz = -d*c1

        d=btable(1,2,0)
        valyx = -d*c1

        d=btable(0,3,0)
        valyy = -d*c1

        d=btable(0,2,1)
        valyz = -d*c1

        d=btable(1,1,1)
        valzx = -d*c1

        d=btable(0,2,1)
        valzy = -d*c1

        d=btable(0,1,2)
        valzz = -d*c1

        d=rtable(1,0,0)
        valyx = valyx+d*2

        d=rtable(0,1,0)
        valyy = valyy+d*2

        d=rtable(0,0,1)
        valyz = valyz+d*2
        endif

        endif

        if ( k .eq. 3 ) then
        
        d=btable(1,0,1)
        valx = -d*c1

        d=btable(0,1,1)
        valy = -d*c1

        d=btable(0,0,2)
        valz = -d*c1

        d=rtable(0,0,0)
        valz = valz+d*2

        if( ifstrain .eq. 1 ) then
        d=btable(2,0,1)
        valxx = -d*c1

        d=btable(1,1,1)
        valxy = -d*c1

        d=btable(1,0,2)
        valxz = -d*c1

        d=btable(1,1,1)
        valyx = -d*c1

        d=btable(0,2,1)
        valyy = -d*c1

        d=btable(0,1,2)
        valyz = -d*c1

        d=btable(1,0,2)
        valzx = -d*c1

        d=btable(0,1,2)
        valzy = -d*c1

        d=btable(0,0,3)
        valzz = -d*c1

        d=rtable(1,0,0)
        valzx = valzx+d*2

        d=rtable(0,1,0)
        valzy = valzy+d*2

        d=rtable(0,0,1)
        valzz = valzz+d*2
        endif

        endif

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
        if( ifstrain .eq. 1 ) then
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2)
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)

        endif

        enddo
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo
c
c
        if( ifstrain .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine elust3triadirecttarg_one_slow
     $     (rlam,rmu,triangle,sigma_sl,ifself,target,ptfrc,strain)
        implicit real *8 (a-h,o-z)
c
c
c     Direct evaluation of displacement and strain due to constant
c     single layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET
c     due to piecewise-constant single layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3)          SLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        dimension triangle(3,3),sigma_sl(3)
        dimension ptfrc(3),strain(3,3)
        dimension w(20),vert1(3),vert2(3),vert3(3)
        dimension vectout(3),vertout(3), ders(3,3)
c
        c1 = (rlam+rmu)/(rlam+2*rmu)
        c2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
        call tri_for_vect(w,sigma_sl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0
c
        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        call triabtable(2,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = -d*c1

        call triabtable(1,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = -d*c1

        call triabtable(1,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = -d*c1

        call triartable(0,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = valx+d*2


        call triabtable(3,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = valxx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = valxy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        call triabtable(1,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = -d*c1

        call triabtable(0,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = -d*c1

        call triabtable(0,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = -d*c1

        call triartable(0,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = valy+d*2


        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(0,3,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = valyx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = valyy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        call triabtable(1,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = -d*c1

        call triabtable(0,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = -d*c1

        call triabtable(0,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = -d*c1

        call triartable(0,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = valz+d*2


        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(0,0,3,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = valzx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = valzy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = valzz+d*2


        endif

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2)
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)

        enddo
c
c
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo

        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine elust3triadirecttarg
     $     (rlam,rmu,ntri,triangles,sigma_sl,
     1     target,ptfrc,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c     Single layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET not lying
c     on the surface due to piecewise-constant single layer density on
c     collection of triangles.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3,ntri)     array of SLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     target(3)            target location
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
        dimension triangles(3,3,1),sigma_sl(3,1)
        dimension ptfrc0(3),strain0(3,3)
        dimension ptfrc(3),strain(3,3),target(3)


        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        do k=1,ntri

        ifself=0
        call elust3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_sl(1,k),
     1     ifself,target,ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo

        enddo

        return
        end
c
c
c
c
c
        subroutine eltst3triadirecttarg_one
     $     (rlam,rmu,triangle,sigma_dl,trinorm,
     $     ifself,target,ptfrc,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c
c     Direct evaluation of displacement and strain due to constant
c     double layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET 
c     due to piecewise-constant double layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3)          DLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     trianorm(3)          triangle normal
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        dimension triangle(3,3),sigma_dl(3),trinorm(3)
        dimension ptfrc(3),strain(3,3)
        dimension w(20),vert1(3),vert2(3),vert3(3)
        dimension vectout(3),vertout(3), ders(3,3)
        dimension tmatr(3,3),umatr(3,3,3)
        dimension smatr(3,3,3),dmatr(3,3,3,3)
        dimension rtable(0:2,0:2,0:2),btable(0:4,0:4,0:4)
c
        C1 = (rlam+rmu)/(rlam+2*rmu)
        C2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
ccc        call tri_for_vect(w,trinorm,vectout)
ccc        call prin2('inside eltst3triadirecttarg, trinorm=*',vectout,3)

        call tri_for_vect(w,sigma_dl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0

        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3
        do i=1,3
        do j=1,3
        umatr(i,j,k)=0
        enddo
        enddo
        enddo
        
        do m=1,3
        do k=1,3
        do i=1,3
        do j=1,3
        dmatr(i,j,k,m)=0
        enddo
        enddo
        enddo
        enddo


        if( ifstrain .eq. 0 ) maxb=3
        if( ifstrain .eq. 0 ) maxr=1
        if( ifstrain .eq. 1 ) maxb=4
        if( ifstrain .eq. 1 ) maxr=2

        do k=0,maxb
        do i=0,maxb
        do j=0,maxb
        if( i+j+k .le. maxb ) then
        call triabtable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        btable(i,j,k)=d
        endif
        enddo
        enddo
        enddo
        
        do k=0,maxr
        do i=0,maxr
        do j=0,maxr
        if( i+j+k .le. maxr ) then
        call triartable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        rtable(i,j,k)=d
        endif
        enddo
        enddo
        enddo


        
        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(3,0,0)
        valxx = -d*c1

        d=btable(2,1,0)
        valxy = -d*c1

        d=btable(2,0,1)
        valxz = -d*c1

        d=btable(2,1,0)
        valyx = -d*c1

        d=btable(1,2,0)
        valyy = -d*c1

        d=btable(1,1,1)
        valyz = -d*c1

        d=btable(2,0,1)
        valzx = -d*c1

        d=btable(1,1,1)
        valzy = -d*c1

        d=btable(1,0,2)
        valzz = -d*c1

        d=rtable(1,0,0)
        valxx = valxx+d*2

        d=rtable(0,1,0)
        valxy = valxy+d*2

        d=rtable(0,0,1)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        d=btable(2,1,0)
        valxx = -d*c1

        d=btable(1,2,0)
        valxy = -d*c1

        d=btable(1,1,1)
        valxz = -d*c1

        d=btable(1,2,0)
        valyx = -d*c1

        d=btable(0,3,0)
        valyy = -d*c1

        d=btable(0,2,1)
        valyz = -d*c1

        d=btable(1,1,1)
        valzx = -d*c1

        d=btable(0,2,1)
        valzy = -d*c1

        d=btable(0,1,2)
        valzz = -d*c1

        d=rtable(1,0,0)
        valyx = valyx+d*2

        d=rtable(0,1,0)
        valyy = valyy+d*2

        d=rtable(0,0,1)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        d=btable(2,0,1)
        valxx = -d*c1

        d=btable(1,1,1)
        valxy = -d*c1

        d=btable(1,0,2)
        valxz = -d*c1

        d=btable(1,1,1)
        valyx = -d*c1

        d=btable(0,2,1)
        valyy = -d*c1

        d=btable(0,1,2)
        valyz = -d*c1

        d=btable(1,0,2)
        valzx = -d*c1

        d=btable(0,1,2)
        valzy = -d*c1

        d=btable(0,0,3)
        valzz = -d*c1

        d=rtable(1,0,0)
        valzx = valzx+d*2

        d=rtable(0,1,0)
        valzy = valzy+d*2

        d=rtable(0,0,1)
        valzz = valzz+d*2

        endif

        umatr(1,k,1)=valxx
        umatr(1,k,2)=valxy
        umatr(1,k,3)=valxz
        umatr(2,k,1)=valyx
        umatr(2,k,2)=valyy
        umatr(2,k,3)=valyz
        umatr(3,k,1)=valzx
        umatr(3,k,2)=valzy
        umatr(3,k,3)=valzz
c
        enddo
c
c
        if( ifstrain .eq. 1 ) then
c       ... and the derivatives 
c
        do m=1,3
c
        if( m .eq. 1 ) then
        ix=1
        iy=0
        iz=0
        endif
        if( m .eq. 2 ) then
        ix=0
        iy=1
        iz=0
        endif
        if( m .eq. 3 ) then
        ix=0
        iy=0
        iz=1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(3+ix,0+iy,0+iz)
        valxx = -d*c1

        d=btable(2+ix,1+iy,0+iz)
        valxy = -d*c1

        d=btable(2+ix,0+iy,1+iz)
        valxz = -d*c1

        d=btable(2+ix,1+iy,0+iz)
        valyx = -d*c1

        d=btable(1+ix,2+iy,0+iz)
        valyy = -d*c1

        d=btable(1+ix,1+iy,1+iz)
        valyz = -d*c1

        d=btable(2+ix,0+iy,1+iz)
        valzx = -d*c1

        d=btable(1+ix,1+iy,1+iz)
        valzy = -d*c1

        d=btable(1+ix,0+iy,2+iz)
        valzz = -d*c1

        d=rtable(1+ix,0+iy,0+iz)
        valxx = valxx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valxy = valxy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        d=btable(2+ix,1+iy,0+iz)
        valxx = -d*c1

        d=btable(1+ix,2+iy,0+iz)
        valxy = -d*c1

        d=btable(1+ix,1+iy,1+iz)
        valxz = -d*c1
        
        d=btable(1+ix,2+iy,0+iz)
        valyx = -d*c1

        d=btable(0+ix,3+iy,0+iz)
        valyy = -d*c1

        d=btable(0+ix,2+iy,1+iz)
        valyz = -d*c1

        d=btable(1+ix,1+iy,1+iz)
        valzx = -d*c1

        d=btable(0+ix,2+iy,1+iz)
        valzy = -d*c1

        d=btable(0+ix,1+iy,2+iz)
        valzz = -d*c1

        d=rtable(1+ix,0+iy,0+iz)
        valyx = valyx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valyy = valyy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        d=btable(2+ix,0+iy,1+iz)
        valxx = -d*c1

        d=btable(1+ix,1+iy,1+iz)
        valxy = -d*c1

        d=btable(1+ix,0+iy,2+iz)
        valxz = -d*c1

        d=btable(1+ix,1+iy,1+iz)
        valyx = -d*c1

        d=btable(0+ix,2+iy,1+iz)
        valyy = -d*c1

        d=btable(0+ix,1+iy,2+iz)
        valyz = -d*c1

        d=btable(1+ix,0+iy,2+iz)
        valzx = -d*c1

        d=btable(0+ix,1+iy,2+iz)
        valzy = -d*c1

        d=btable(0+ix,0+iy,3+iz)
        valzz = -d*c1

        d=rtable(1+ix,0+iy,0+iz)
        valzx = valzx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valzy = valzy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valzz = valzz+d*2

        endif

        dmatr(1,k,1,m)=valxx
        dmatr(1,k,2,m)=valxy
        dmatr(1,k,3,m)=valxz
        dmatr(2,k,1,m)=valyx
        dmatr(2,k,2,m)=valyy
        dmatr(2,k,3,m)=valyz
        dmatr(3,k,1,m)=valzx
        dmatr(3,k,2,m)=valzy
        dmatr(3,k,3,m)=valzz
c
        enddo
        enddo
c
ccc        call prin2('umatr=*',umatr,3*3*3)
c
        endif
c
c
        do i=1,3
        do j=1,3
        tmatr(i,j)=0
        enddo
        enddo
c
        do i=1,3
        tmatr(i,3)=rlam*(umatr(1,i,1)+umatr(2,i,2)+umatr(3,i,3)) 
        enddo
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+rmu*(umatr(j,i,3)+umatr(3,i,j))
        enddo
        enddo
c
ccc        call prin2('tmatr=*',tmatr,3*3)
c
        if( ifstrain .eq. 1 ) then
c
        do m=1,3
c
        do i=1,3
        do j=1,3
        smatr(i,j,m)=0
        enddo
        enddo
c
        do i=1,3
        smatr(i,3,m)=
     $     rlam*(dmatr(1,i,1,m)+dmatr(2,i,2,m)+dmatr(3,i,3,m)) 
        enddo
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+rmu*(dmatr(j,i,3,m)+dmatr(3,i,j,m))
        enddo
        enddo
c
        enddo
c
ccc        call prin2('smatr=*',tmatr,3*3*3)
c
        endif
c
        do k=1,3

        valx=tmatr(1,k)
        valy=tmatr(2,k)
        valz=tmatr(3,k)

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
        if( ifstrain .eq. 1 ) then
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxx=smatr(1,k,1)
        valyx=smatr(2,k,1)
        valzx=smatr(3,k,1)
        valxy=smatr(1,k,2)
        valyy=smatr(2,k,2)
        valzy=smatr(3,k,2)
        valxz=smatr(1,k,3)
        valyz=smatr(2,k,3)
        valzz=smatr(3,k,3)

        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2) 
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)
c
        endif
c
        enddo
c
c
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo
c
        if( ifstrain .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine eltst3triadirecttarg_one_slow
     $     (rlam,rmu,triangle,sigma_dl,trinorm,
     $     ifself,target,ptfrc,strain)
        implicit real *8 (a-h,o-z)
c
c
c     Direct evaluation of displacement and strain due to constant
c     double layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET 
c     due to piecewise-constant double layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3)          DLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     trianorm(3)          triangle normal
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        dimension triangle(3,3),sigma_dl(3),trinorm(3)
        dimension ptfrc(3),strain(3,3)
        dimension w(20),vert1(3),vert2(3),vert3(3)
        dimension vectout(3),vertout(3), ders(3,3)
        dimension tmatr(3,3),umatr(3,3,3)
        dimension smatr(3,3,3),dmatr(3,3,3,3)
c
        c1 = (rlam+rmu)/(rlam+2*rmu)
        c2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
ccc        call tri_for_vect(w,trinorm,vectout)
ccc        call prin2('inside eltst3triadirecttarg, trinorm=*',vectout,3)

        call tri_for_vect(w,sigma_dl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0

        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3
        do i=1,3
        do j=1,3
        umatr(i,j,k)=0
        enddo
        enddo
        enddo
        
        do m=1,3
        do k=1,3
        do i=1,3
        do j=1,3
        dmatr(i,j,k,m)=0
        enddo
        enddo
        enddo
        enddo
        
        do k=1,3

        if ( k .eq. 1 ) then
        
        call triabtable(3,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = valxx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = valxy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(0,3,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = valyx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = valyy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(0,0,3,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = valzx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = valzy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = valzz+d*2

        endif

        umatr(1,k,1)=valxx
        umatr(1,k,2)=valxy
        umatr(1,k,3)=valxz
        umatr(2,k,1)=valyx
        umatr(2,k,2)=valyy
        umatr(2,k,3)=valyz
        umatr(3,k,1)=valzx
        umatr(3,k,2)=valzy
        umatr(3,k,3)=valzz
c
        enddo
c
c
c       ... and the derivatives 
c
        do m=1,3
c
        if( m .eq. 1 ) then
        ix=1
        iy=0
        iz=0
        endif
        if( m .eq. 2 ) then
        ix=0
        iy=1
        iz=0
        endif
        if( m .eq. 3 ) then
        ix=0
        iy=0
        iz=1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        call triabtable(3+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(2+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(2+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(2+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(1+ix,2+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(2+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(1+ix,0+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = valxx+d*2

        call triartable(0+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = valxy+d*2

        call triartable(0+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        call triabtable(2+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(1+ix,2+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(1+ix,2+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(0+ix,3+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(0+ix,2+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(0+ix,2+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(0+ix,1+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = valyx+d*2

        call triartable(0+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = valyy+d*2

        call triartable(0+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        call triabtable(2+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*c1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*c1

        call triabtable(1+ix,0+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*c1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*c1

        call triabtable(0+ix,2+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*c1

        call triabtable(0+ix,1+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*c1

        call triabtable(1+ix,0+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*c1

        call triabtable(0+ix,1+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*c1

        call triabtable(0+ix,0+iy,3+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*c1

        call triartable(1+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = valzx+d*2

        call triartable(0+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = valzy+d*2

        call triartable(0+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = valzz+d*2

        endif

        dmatr(1,k,1,m)=valxx
        dmatr(1,k,2,m)=valxy
        dmatr(1,k,3,m)=valxz
        dmatr(2,k,1,m)=valyx
        dmatr(2,k,2,m)=valyy
        dmatr(2,k,3,m)=valyz
        dmatr(3,k,1,m)=valzx
        dmatr(3,k,2,m)=valzy
        dmatr(3,k,3,m)=valzz
c
        enddo
        enddo
c
c
ccc        call prin2('umatr=*',umatr,3*3*3)
c
c
        do i=1,3
        do j=1,3
        tmatr(i,j)=0
        enddo
        enddo
c
        do i=1,3
        tmatr(i,3)=rlam*(umatr(1,i,1)+umatr(2,i,2)+umatr(3,i,3)) 
        enddo
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+rmu*(umatr(j,i,3)+umatr(3,i,j))
        enddo
        enddo
c
ccc        call prin2('tmatr=*',tmatr,3*3)
c
        do m=1,3
c
        do i=1,3
        do j=1,3
        smatr(i,j,m)=0
        enddo
        enddo
c
        do i=1,3
        smatr(i,3,m)=
     $     rlam*(dmatr(1,i,1,m)+dmatr(2,i,2,m)+dmatr(3,i,3,m)) 
        enddo
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+rmu*(dmatr(j,i,3,m)+dmatr(3,i,j,m))
        enddo
        enddo
c
        enddo
c
ccc        call prin2('smatr=*',tmatr,3*3*3)
c
c
        do k=1,3

        valx=tmatr(1,k)
        valy=tmatr(2,k)
        valz=tmatr(3,k)

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxx=smatr(1,k,1)
        valyx=smatr(2,k,1)
        valzx=smatr(3,k,1)
        valxy=smatr(1,k,2)
        valyy=smatr(2,k,2)
        valzy=smatr(3,k,2)
        valxz=smatr(1,k,3)
        valyz=smatr(2,k,3)
        valzz=smatr(3,k,3)

        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2) 
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)

        enddo
c
c
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo

        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine eltst3triadirecttarg
     $     (rlam,rmu,ntri,triangles,sigma_dl,trinorm,
     1     target,ptfrc,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c     Double layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET not lying
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3,ntri)     array of DLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     trianorm(3,ntri)     array of triangle normals
c     target(3)            target location
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
        dimension triangles(3,3,1),sigma_dl(3,1),trinorm(3,1)
        dimension ptfrc0(3),strain0(3,3)
        dimension ptfrc(3),strain(3,3),target(3)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        do k=1,ntri

        ifself=0
        call eltst3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_dl(1,k),trinorm(1,k),
     1     ifself,target,ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
        subroutine elust3triadirectself
     $     (rlam,rmu,ipatch,ntri,triangles,sigma_sl,
     1     zparts,ptfrc,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c     Single layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at centroid zparts(*,ipatch) 
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles, numbered jpatch = 1,...,ntri. 
c
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3,ntri)     array of SLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     zparts(3,1)          array of triangle centroids
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at centroid zparts(*,ipatch)
c     strain(3,3)         strain at centroid zparts(*,ipatch)
c
c
        dimension triangles(3,3,1),sigma_sl(3,1),trinorm(3,1)
        dimension ptfrc0(3),strain0(3,3)
        dimension ptfrc(3),strain(3,3),zparts(3,1)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        do k=1,ntri

        if( k .eq. ipatch ) ifself=1
        if( k .ne. ipatch ) ifself=0

        call elust3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_sl(1,k),
     1     ifself,zparts(1,ipatch),ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
        subroutine eltst3triadirectself
     $     (rlam,rmu,ipatch,ntri,triangles,sigma_dl,trinorm,
     1     zparts,ptfrc,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c     Double layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at centroid zparts(*,ipatch) 
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles, numbered jpatch = 1,...,ntri. 
c
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3,ntri)     array of DLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     trianorm(3,ntri)     array of triangle normals
c     zparts(3,1)          array of triangle centroids
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at centroid zparts(*,ipatch)
c     strain(3,3)         strain at centroid zparts(*,ipatch)
c
c
        dimension triangles(3,3,1),sigma_dl(3,1),trinorm(3,1)
        dimension ptfrc0(3),strain0(3,3)
        dimension ptfrc(3),strain(3,3),zparts(3,1)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo
        
        do k=1,ntri

        if( k .eq. ipatch ) ifself=1
        if( k .ne. ipatch ) ifself=0

        call eltst3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_dl(1,k),trinorm(1,k),
     1     ifself,zparts(1,ipatch),ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
