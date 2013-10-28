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
c        potentials in free space in R^3.  Half space elastostatic
c        Green's function (Mindlin's solution).        
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      elhfmm3dtriatargiter - evaluates the elastostatic potential ON OR OFF
c         SURFACE due to a collection of flat triangles with constant
c         single and/or double layer densities using the Fast Multipole
c         Method (Mindlin's solution).
c
c      elh3dtriadirecttarg - evaluates the elastostatic potential ON OR
c         OFF SURFACE due to a collection of flat triangles with
c         constant single and/or double layer densities using the direct
c         O(N^2) algorithm (Mindlin's solution).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine elhfmm3dtriatargiter
     $     (ier,iprec,RLAM,RMU,TRIANGLE,TRINORM,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg,
     $     icomp_type,wsave,lwsave,lused)
C
c     This is the principal subroutine for evaluating elastostatic layer
c     potentials on (flat) triangulated surfaces.  It permits the
c     evaluation of a single layer potential with piecewise constant
c     density defined by the (force) vector sigma_sl and a dipole layer
c     potential with piecewise constant density vector sigma_dl and
c     dipole orientation defined by trinorm.
c
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c     Half space elastostatic Green's function (Mindlin's solution).
c
c
c     NOTE #1: In the present version, dipole vectors SIGMA_DV are SET to
c     be EQUAL to the triangle normals.
c
c     NOTE #2: Precomputation of local part for half space image FMM is
c     meaningful for targets only. Currently, self evaluation should be
c     done with a separate call, by passing sources as targets. If both
c     on and off surface evaluation is required, sources and targets can
c     be merged into one joined target list.
c
c
c 
C     INPUT:
C
C     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     nparts = number of sources
C     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
C     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
C     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
C     target(3,ntargs) = evaluation target points
c     ifptfrctarg - target displacement computation flag
c     ifstraintarg - target strain computation flag
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
c     icomp_type: FMM precomputation flag
c        icomp_type = 1, the FMM is carried out without local storage.
c        icomp_type = 2, the FMM is carried out AND it stores
c               local  interactions in wsave_direct storage array.
c        icomp_type = 3, the FMM is carried USING the local 
c               interactions stored in wsave_direct. It assumes a previous
c               call with icomp_type=2 was successfully completed.       
c               It requires the same FMM tree, so must be called with the same
c               parameters as the previous call to the FMM.
c        icomp_type = 4 estimates the amount of workspace needed 
c               for the storage of local interactions and sets lused parameter.
c               It doesn't actually do the FMM calculation, except for tree
c               generation.
c
C     OUTPUT:
C
C     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
c
c     ier - error return code
c             ier=0     =>  normal execution
c             ier=4     =>  cannot allocate tree workspace
c             ier=8     =>  cannot allocate bulk FMM  workspace
c             ier=16    =>  cannot allocate mpole expansion workspace in FMM
c
c
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,nparts),trinorm(3,nparts)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)

ccc        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 ptfrc(3,1),strain(3,3,1)

        real *8 target(3,ntargs)
        real *8 ptfrcTARG(3,ntargs),strainTARG(3,3,ntargs)
        integer nparts,ntargs
c
        real *8, allocatable :: sourceim(:,:)
        real *8, allocatable :: targetim(:,:)
c
        real *8, allocatable :: source_interp(:,:,:)
        real *8, allocatable :: sigma_sl_interp(:,:,:)
        real *8, allocatable :: sigma_dl_interp(:,:,:)
        real *8, allocatable :: sigma_dv_interp(:,:,:)
c
        real *8, allocatable :: ptfrc1(:,:)
        real *8, allocatable :: ptfrctarg1(:,:)
        real *8, allocatable :: strain1(:,:,:)
        real *8, allocatable :: straintarg1(:,:,:)
c
        dimension vert1(3),vert2(3),vert3(3),w(1000)
        dimension vertout(3)
        dimension vert01(3),vert02(3),vert03(3)
        dimension rnodes(2,1000)
        dimension weights(1000)
c
        dimension ptfrc0(3),strain0(3,3)
c

        real *8 wsave(1)
c
c
c       ... disable self evaluation routines 
c       for both displacements and strains
c
        ifptfrc=0
        ifstrain=0
c
c
c
        if( ifptfrc .eq. 1 ) then
        do i=1,nparts
        ptfrc(1,i) = 0
        ptfrc(2,i) = 0
        ptfrc(3,i) = 0
        enddo
        endif
c
        if( ifstrain .eq. 1 ) then
        do i=1,nparts
        do j=1,3
        do k=1,3
        strain(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c

        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
        ptfrctarg(1,i) = 0
        ptfrctarg(2,i) = 0
        ptfrctarg(3,i) = 0
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
        do j=1,3
        do k=1,3
        straintarg(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c
c
        ifevalfar=1
        ifevalloc=1
c
        if( icomp_type .eq. 4 .or. icomp_type .eq. 5 ) then
        ifevalfar=0
        ifevalloc=1
        endif
c
        if( icomp_type .eq. 0 ) then
        ifevalfar=0
        ifevalloc=0
        endif
c
c
        if( ifevalfar .eq. 1 ) then
        call elhfmm3dtriatargiter_main
     $     (ier,iprec,RLAM,RMU,TRIANGLE,TRINORM,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg,
     $     icomp_type,wsave,lwsave,lused)
        endif
c
c       ... disable local evaluation for itype = 1,2
c       return
c
c
c     reconstruct FMM data structure and account for  all local 
c     interactions using quadrature routines for piecewise
c     constant densities.
c
c
        if( ifevalloc .eq. 1 ) then

        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
        allocate( targetim(3,ntargs) )
        do i=1,ntargs
        targetim(1,i)= target(1,i)
        targetim(2,i)= target(2,i)
        targetim(3,i)=-target(3,i)
        enddo
        endif
c
        call elhfmm3dtriatargiter0
     $     (ier,iprec,RLAM,RMU,TRIANGLE,TRINORM,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,TRINORM,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,targetim,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg,
     $     icomp_type,wsave,lwsave,lused)

        endif
c
        return
        end
c
c
c
c
c
C*********************************
      subroutine elhfmm3dtriatargiter_main
     $     (ier,iprec,RLAM,RMU,TRIANGLE,TRINORM,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg,
     $     icomp_type,wsave,lwsave,lused)
C
c     This is the principal subroutine for evaluating elastostatic layer
c     potentials on (flat) triangulated surfaces.  It permits the
c     evaluation of a single layer potential with piecewise constant
c     density defined by the (force) vector sigma_sl and a dipole layer
c     potential with piecewise constant density vector sigma_dl and
c     dipole orientation defined by trinorm.
c
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c     Half space elastostatic Green's function (Mindlin's solution).
c
c
c     NOTE: In the present version, dipole vectors SIGMA_DV are SET to
c     be EQUAL to the triangle normals.
c
c 
C     INPUT:
C
C     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     nparts = number of sources
C     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
C     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
C     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
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
C
c
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,nparts),trinorm(3,nparts)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrcTARG(3,ntargs),strainTARG(3,3,ntargs)
        integer nparts,ntargs
c
        real *8, allocatable :: source_interp(:,:,:)
        real *8, allocatable :: sigma_sl_interp(:,:,:)
        real *8, allocatable :: sigma_dl_interp(:,:,:)
        real *8, allocatable :: sigma_dv_interp(:,:,:)
c
        real *8, allocatable :: ptfrc1(:,:)
        real *8, allocatable :: ptfrctarg1(:,:)
        real *8, allocatable :: strain1(:,:,:)
        real *8, allocatable :: straintarg1(:,:,:)
c
        dimension vert1(3),vert2(3),vert3(3),w(1000)
        dimension vertout(3)
        dimension vert01(3),vert02(3),vert03(3)
        dimension rnodes(2,1000)
        dimension weights(1000)
c
        dimension ptfrc0(3),strain0(3,3)
c
        real *8 wsave(1)
c
        if( ifptfrc .eq. 1 ) then
        do i=1,nparts
        ptfrc(1,i) = 0
        ptfrc(2,i) = 0
        ptfrc(3,i) = 0
        enddo
        endif
c
        if( ifstrain .eq. 1 ) then
        do i=1,nparts
        do j=1,3
        do k=1,3
        strain(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c

        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
        ptfrctarg(1,i) = 0
        ptfrctarg(2,i) = 0
        ptfrctarg(3,i) = 0
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
        do j=1,3
        do k=1,3
        straintarg(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c
c
c***********************************************************************
c     PART 1 === direct arrival
c***********************************************************************
c
c***********************************************************************
c     PART 2,3,4 === image contribution
c***********************************************************************
c
c
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
        allocate( ptfrc1(3,nparts) )
        allocate( strain1(3,3,nparts) )
        else
        allocate( ptfrc1(3,1) )
        allocate( strain1(3,3,1) )
        endif
c
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
        allocate( ptfrctarg1(3,ntargs) )
        allocate( straintarg1(3,3,ntargs) )
        else
        allocate( ptfrctarg1(3,1) )
        allocate( straintarg1(3,3,1) )
        endif
c
c
        if( ifptfrc .eq. 1 ) then
        do i=1,nparts
        ptfrc1(1,i) = 0
        ptfrc1(2,i) = 0
        ptfrc1(3,i) = 0
        enddo
        endif
c
        if( ifstrain .eq. 1 ) then
        do i=1,nparts
        do j=1,3
        do k=1,3
        strain1(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c

        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
        ptfrctarg1(1,i) = 0
        ptfrctarg1(2,i) = 0
        ptfrctarg1(3,i) = 0
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
        do j=1,3
        do k=1,3
        straintarg1(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c
c
c       itype = 1: use particle fmm (oversampled)
c       itype = 2: adaptive integration
c       itype = 3: use triangle fmm
c
c
        itype = 3

c
        if( itype .eq. 3 ) then
c
c       ... process sources and targets separately
c
        ifptfrc0=0
        ifstrain0=0

        call elhfmm3dtriatarg_image
     $     (ier,iprec,RLAM,RMU,triangle,trinorm,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,
     $     ifdouble,SIGMA_DL,TRINORM,
     $     ifptfrc0,ptfrc1,ifstrain0,strain1,
     $     nparts,source,ifptfrc,PTFRC1,
     $     ifstrain,STRAIN1)
c
        call elhfmm3dtriatarg_image
     $     (ier,iprec,RLAM,RMU,triangle,trinorm,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,
     $     ifdouble,SIGMA_DL,TRINORM,
     $     ifptfrc0,ptfrc1,ifstrain0,strain1,
     $     ntargs,target,ifptfrctarg,PTFRCtarg1,
     $     ifstraintarg,STRAINtarg1)
c
        goto 8000
        endif
c
c
c       
 8000   continue
c
c
        if( ifptfrc .eq. 1 ) then
        do i=1,nparts
        ptfrc(1,i) = ptfrc(1,i) + ptfrc1(1,i)
        ptfrc(2,i) = ptfrc(2,i) + ptfrc1(2,i)
        ptfrc(3,i) = ptfrc(3,i) + ptfrc1(3,i)
        enddo
        endif
c
        if( ifstrain .eq. 1 ) then
        do i=1,nparts
        do j=1,3
        do k=1,3
        strain(j,k,i) = strain(j,k,i) + strain1(j,k,i)
        enddo
        enddo
        enddo
        endif
c

        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
        ptfrctarg(1,i) = ptfrctarg(1,i) + ptfrctarg1(1,i)
        ptfrctarg(2,i) = ptfrctarg(2,i) + ptfrctarg1(2,i)
        ptfrctarg(3,i) = ptfrctarg(3,i) + ptfrctarg1(3,i)
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
        do j=1,3
        do k=1,3
        straintarg(j,k,i) = straintarg(j,k,i) + straintarg1(j,k,i)
        enddo
        enddo
        enddo
        endif
c
        return
        end
c
c
c
c
c
      SUBROUTINE elhfmm3dtriatarg_image
     $     (ier,iprec,RLAM,RMU,triangle,trinorm,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)
C
C     FMM calculation subroutine for elastostatic N-body problem.
C
C     Half space elastostatic Green's function (Mindlin's solution).
C
C     INPUT:
C
C     RLAM, RMU = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     NPARTS = number of sources
C     SOURCE(3,NPARTS) = source locations
c     ifsingle = single layer computation flag  
C     SIGMA_SL(3,NPARTS) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
C     SIGMA_DL(3,NPARTS) = vector strength of nth dipole (double layer)
C     SIGMA_DV(3,NPARTS) = dipole orientation vectors (double layer)
C     TARGET(3,NTARGS) = evaluation target points
C
c     iprec:  FMM precision flag
c
C     OUTPUT:
C
C     PTFRC(3,NPARTS) = displacement 
c     STRAIN(3,3,NPARTS) = strain
C     PTFRCtarg(3,NTARGS) = displacement at target locations
c     STRAINtarg(3,3,NTARGS) = strain at target locations
C
c
        IMPLICIT REAL *8 (A-H,O-Z)
        real *8 triangle(3,3,nparts),trinorm(3,nparts)
        REAL *8 SOURCE(3,NPARTS)
        REAL *8 SIGMA_SL(3,NPARTS)
        REAL *8 SIGMA_DL(3,NPARTS),SIGMA_DV(3,NPARTS)
        REAL *8 PTFRC(3,NPARTS),STRAIN(3,3,NPARTS)
        REAL *8 TARGET(3,NTARGS)
        REAL *8 PTFRCtarg(3,NTARGS),STRAINtarg(3,3,NTARGS)
        INTEGER NPARTS,NTARGS
        real *8, allocatable :: w(:)
c       
        real *8, allocatable :: sourceim(:,:)
        real *8, allocatable :: targetim(:,:)
c
        real *8, allocatable :: ptfrc1(:,:)
        real *8, allocatable :: ptfrctarg1(:,:)
        real *8, allocatable :: strain1(:,:,:)
        real *8, allocatable :: straintarg1(:,:,:)
c
        real *8, allocatable :: sigma_sl_im(:,:)
        real *8, allocatable :: sigma_dl_im(:,:)
        real *8, allocatable :: sigma_dv_im(:,:)
        real *8, allocatable :: triangleim(:,:,:)
        real *8, allocatable :: trinormim(:,:)
c
        complex *16, allocatable :: charge(:)
        complex *16, allocatable :: pottarg(:,:)
        complex *16, allocatable :: fldtarg(:,:,:)
c
        real *8 rlame(2)
c
c
        if( ifptfrc .eq. 1 ) then
        do i=1,nparts
        ptfrc(1,i) = 0
        ptfrc(2,i) = 0
        ptfrc(3,i) = 0
        enddo
        endif
c
        if( ifstrain .eq. 1 ) then
        do i=1,nparts
        do j=1,3
        do k=1,3
        strain(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c
c
        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
        ptfrctarg(1,i) = 0
        ptfrctarg(2,i) = 0
        ptfrctarg(3,i) = 0
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
        do j=1,3
        do k=1,3
        straintarg(j,k,i) = 0
        enddo
        enddo
        enddo
        endif
c
c
c***********************************************************************
c     PART 1 === direct arrival
c***********************************************************************
c
        ifptfrc0=0
        ifstrain0=0
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stresslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c
c
        allocate( sourceim(3,nparts) )
        allocate( trinormim(3,nparts) )
        allocate( triangleim(3,3,nparts) )
        do i=1,nparts
        sourceim(1,i)= source(1,i)
        sourceim(2,i)= source(2,i)
        sourceim(3,i)=-source(3,i)
        enddo
        do i=1,nparts
        triangleim(1,1,i)= triangle(1,1,i)
        triangleim(2,1,i)= triangle(2,1,i)
        triangleim(3,1,i)= -triangle(3,1,i)
        triangleim(1,2,i)= triangle(1,2,i)
        triangleim(2,2,i)= triangle(2,2,i)
        triangleim(3,2,i)= -triangle(3,2,i)
        triangleim(1,3,i)= triangle(1,3,i)
        triangleim(2,3,i)= triangle(2,3,i)
        triangleim(3,3,i)= -triangle(3,3,i)
c
        trinormim(1,i)= trinorm(1,i)
        trinormim(2,i)= trinorm(2,i)
        trinormim(3,i)=-trinorm(3,i)
        enddo
c
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
        allocate( targetim(3,ntargs) )
        do i=1,ntargs
        targetim(1,i)= target(1,i)
        targetim(2,i)= target(2,i)
        targetim(3,i)=-target(3,i)
        enddo
        allocate( ptfrctarg1(3,ntargs) )
        allocate( straintarg1(3,3,ntargs) )
        endif
c
        ifptfrc0=0
        ifstrain0=0

        if( ifsingle .eq. 1 ) then

        allocate( sigma_sl_im(3,nparts) )
        do i=1,nparts
        sigma_sl_im(1,i)=-sigma_sl(1,i)
        sigma_sl_im(2,i)=-sigma_sl(2,i)
        sigma_sl_im(3,i)=-sigma_sl(3,i)
        enddo

        endif
c

        if( ifdouble .eq. 1 ) then
c
c       we need to flip the sign here, so that we could use one FMM call 
c       if both single and double layers are present
c        
        allocate( sigma_dl_im(3,nparts) )
        do i=1,nparts
        sigma_dl_im(1,i)=-sigma_dl(1,i)
        sigma_dl_im(2,i)=-sigma_dl(2,i)
        sigma_dl_im(3,i)=-sigma_dl(3,i)
        enddo
c
        allocate( sigma_dv_im(3,nparts) )
        do i=1,nparts
        sigma_dv_im(1,i)=-sigma_dv(1,i)
        sigma_dv_im(2,i)=-sigma_dv(2,i)
        sigma_dv_im(3,i)=sigma_dv(3,i)
        enddo
c
        else
c
c       fmm needs a memory reference to proceed, allocate a small array
c
        allocate( sigma_dl_im(3,1) )
c
        endif
c
        if_mindlin_only = 0
c       
        if( if_mindlin_only .eq. 1 ) goto 5000
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c
        RLAMIM = RLAM+4*RMU
        RMUIM = -RMU
        call prinf('===inside elxfmm3dtriatarg===*',i,0)
        call elxfmm3dtriatarg
     $     (ier,iprec,RLAMIM,RMUIM,triangle,trinorm,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL_IM,
     $     ifptfrc0,ptfrc,ifstrain0,strain,
     $     ntargs,targetim,ifptfrctarg,PTFRCtarg1,
     $     ifstraintarg,STRAINtarg1)
c
c     scaling is by 1/(2*rmuim) and it should be by 1/(2*rmu).
c     This is fixed by flipping the sign (subtracting) in contributions
c     to displacement and strain here.
c
        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
        ptfrctarg(1,i) = ptfrctarg(1,i) - ptfrctarg1(1,i)
        ptfrctarg(2,i) = ptfrctarg(2,i) - ptfrctarg1(2,i)
        ptfrctarg(3,i) = ptfrctarg(3,i) - ptfrctarg1(3,i) 
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
        do j=1,3
        do k=1,3
        straintarg(j,k,i) = straintarg(j,k,i) - straintarg1(j,k,i)
        enddo
        enddo
        enddo
        endif
c
        endif
c
ccc        return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source, parts B and C
c***********************************************************************
c
 5000   continue
c
        rlame(1) = rlam
        rlame(2) = rmu
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
        allocate( pottarg(3,ntargs) )
        allocate( fldtarg(3,3,ntargs) )
c
        call prinf('===inside lfmm3dmindlintriatarg===*',i,0)
        ifpottarg=1
        iffldtarg=1
c
        call lfmm3dmindlintriatarg(ier,iprec,
     $     nparts,triangleim,trinormim,sourceim,
     $     rlame,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ntargs,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
        if( ifptfrctarg .eq. 1 ) then
        do icomp = 1,3
        do i=1,ntargs
        ptfrctarg(icomp,i)=ptfrctarg(icomp,i)+pottarg(icomp,i)/rmu
        enddo
        enddo
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do icomp = 1,3
        do i=1,ntargs
        straintarg(1,icomp,i)=
     $     straintarg(1,icomp,i)-fldtarg(1,icomp,i)/(2*rmu)
        straintarg(2,icomp,i)=
     $     straintarg(2,icomp,i)-fldtarg(2,icomp,i)/(2*rmu)
        straintarg(3,icomp,i)=
     $     straintarg(3,icomp,i)-fldtarg(3,icomp,i)/(2*rmu)
        straintarg(icomp,1,i)=
     $     straintarg(icomp,1,i)-fldtarg(1,icomp,i)/(2*rmu)
        straintarg(icomp,2,i)=
     $     straintarg(icomp,2,i)-fldtarg(2,icomp,i)/(2*rmu)
        straintarg(icomp,3,i)=
     $     straintarg(icomp,3,i)-fldtarg(3,icomp,i)/(2*rmu)
        enddo
        enddo
        endif
c
        deallocate( pottarg )
        deallocate( fldtarg )
        endif
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine elxfmm3dtria
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain)
C
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
C     INPUT:
C
C     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     nparts = number of sources
C     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
C     sigma_sl(3,nparts) = vector strength of single layer source 
c     ifdouble = double layer computation flag  
C     sigma_dl(3,nparts) = vector strength of double layer source
C
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
C     OUTPUT:
C
C     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
C
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
        call elxfmm3dtriatarg
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
C*********************************
      subroutine elxfmm3dtriatarg
     $     (ier,iprec,RLAM,RMU,TRIANGLE,TRINORM,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)
C
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
C     INPUT:
C
C     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     nparts = number of sources
C     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
C     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
C     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
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
C
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrcTARG(3,ntargs),strainTARG(3,3,ntargs)
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
        call elxfmm3dtriatargmain_fast
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,trinorm,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg,
     $     w(icharge),w(idipstr),w(idipvec),
     $     w(ipot),w(ifld),w(ihess),w(ihessmatr),
     $     w(ipottarg),w(ifldtarg),w(ihesstarg),w(ihessmatrtarg))
c
        return
        end
c
c
c
c
c
C*********************************
      subroutine elxfmm3dtriatargmain_fast
     $     (ier,iprec,rlam,rmu,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg,
     $     charge,dipstr,dipvec,pot,fld,hess,hessmatr,
     $     pottarg,fldtarg,hesstarg,hessmatrtarg)
C
C     FMM calculation subroutine for elastostatic N-body problem
c
c     4 Laplace FMM calls, uses linear charge and dipole densities
C
c
C     INPUT:
C
C     rlam, rmu = Lame parameters
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
C     nparts = number of sources
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
C     OUTPUT:
C
C     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
C     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
C
c
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,nparts),trinorm(3,nparts)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrctarg(3,ntargs),strainTARG(3,3,ntargs)
        integer nparts,ntargs
c       
        complex *16 charge(3,1)
        complex *16 dipstr(3,1)
        REAL *8 dipvec(3,3,1)
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
           call lfmm3dtrilhesstargiter(ier,iprec,
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
     $       sigma_sl(3,k)*triangle(3,m,k))*(c1)/(2*rmu)
          ifcharge = 1
          endif
          if( ifdouble .eq. 1 ) then
          charge(m,k) = charge(m,k) + 
     $        (sigma_dl(1,k)*sigma_dv(1,k) + 
     1         sigma_dl(2,k)*sigma_dv(2,k) + 
     2         sigma_dl(3,k)*sigma_dv(3,k))*(c2+2)
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
        call lfmm3dtrilhesstargiter(ier,iprec,
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
              PTFRC(1,k)=PTFRC(1,k)
              PTFRC(2,k)=PTFRC(2,k)
              PTFRC(3,k)=PTFRC(3,k)
           endif
c       
           if( ifstrain .eq. 1 ) then
           do i=1,3
           hessmatr(i,3,k)=-hessmatr(i,3,k)
           enddo
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
        hessmatrtarg(i,3,k)=-hessmatrtarg(i,3,k)
        enddo
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
        subroutine elhfmm3dtriatargiter0(ier,iprec,RLAM,RMU,
     $     triangle,trinorm,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,
     $     ifptfrctarg,ptfrctarg,ifstraintarg,straintarg,
     $     icomp_type,wsave,lwsave,lused)
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
        real *8 wsave(1)
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
        bscale = 1/1.5
        bscale = 2
        if( iprec .eq. -2 ) nbox=8*bscale
        if( iprec .eq. -1 ) nbox=15*bscale
        if( iprec .eq. 0 ) nbox=30*bscale
        if( iprec .eq. 1 ) nbox=60*bscale
        if( iprec .eq. 2 ) nbox=120*bscale
        if( iprec .eq. 3 ) nbox=240*bscale
        if( iprec .eq. 4 ) nbox=480*bscale
        if( iprec .eq. 5 ) nbox=700*bscale
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
        call elhfmm3dtriatargiter0_evalloc(ier,iprec,RLAM,RMU,
     $     w(itrianglesort),w(itrinormsort),
     $     nsource,w(isourcesort),
     $     ifsingle,w(isigma_slsort),
     $     ifdouble,w(isigma_dlsort),w(isigma_dvsort),
     $     ifptfrc,w(iptfrc),ifstrain,w(istrain),
     $     ntarget,w(itargetsort),
     $     ifptfrctarg,w(iptfrctarg),ifstraintarg,w(istraintarg),
     $     nboxes,laddr,nlev,wlists(iwlists),lwlists,
     $     icomp_type,wsave,lwsave,lused)
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
        subroutine elhfmm3dtriatargiter0_evalloc(ier,iprec,RLAM,RMU,
     $     triangle,trinorm,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,
     $     ifptfrctarg,ptfrctarg,ifstraintarg,straintarg,
     $     nboxes,laddr,nlev,wlists,lwlists,
     $     icomp_type,wsave,lwsave,lused)
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
        integer, allocatable :: iwboxes(:,:)
        integer, allocatable :: iwpairs(:,:)
c
        real *8 wsave(1)
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

        ifprint = 0

ccc        if( ifevalloc .eq. 0 ) goto 9000
c
        call prinf('=== STEP 8a (count) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8a, count direct interactions 
c
c
        allocate( iwboxes(4,nboxes) )
c
c
c       ... count number of interaction pairs
c
        npairs_count = 0
c
        do 6100 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
c        if (nkids .eq. 0) then
c            npts=box(15)
c            if (ifprint .ge. 2) then
c               call prinf('npts=*',npts,1)
c               call prinf('isource=*',isource(box(14)),box(15))
c            endif
c        endif
c
c       ... prune all sourceless boxes
c       adjust processing for lists 3 and 4
c
ccc        if (nkids .eq. 0 .and. box(15) .ne. 0 ) then
        if ( box(15) .ne. 0 ) then
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist1,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist1)
c
        nlist = nlist1
c
c
c       ... evaluate self interactions
c
         if( nkids .eq. 0 ) then
         list(nlist+1) = ibox
         nlist = nlist + 1 
         endif
c    
c
c       ... evaluate interactions with list #3 
c
        itype=3
        call d3tgetl(ier,ibox,itype,list(nlist+1),nlist3,wlists)
        if (ifprint .ge. 2) call prinf('list3=*',list,nlist3)
        nlist = nlist + nlist3
c
c    
c       ... evaluate interactions with list #4
c
        itype=4
        call d3tgetl(ier,ibox,itype,list(nlist+1),nlist4,wlists)
        if (ifprint .ge. 2) call prinf('list4=*',list,nlist4)
        nlist = nlist + nlist4
c
c
c       ... for all pairs in list #1, #3, and #4
c       evaluate the potentials and fields directly
c
        npairs_count = npairs_count + nlist
c
        endif
c
 6100   continue
c
c
        allocate( iwpairs(4,npairs_count) )
c
        npairs = 0
        ntotal = 0
        iptr = 1
        lptr = 1
        next_pair  = 1 
        next_chunk = 1
c
c
c
        do 6200 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
c        if (nkids .eq. 0) then
c            npts=box(15)
c            if (ifprint .ge. 2) then
c               call prinf('npts=*',npts,1)
c               call prinf('isource=*',isource(box(14)),box(15))
c            endif
c        endif
c
        lpairs=0
        lstore=0
c
c       ... prune all sourceless boxes
c       adjust processing for lists 3 and 4
c
ccc        if (nkids .eq. 0 .and. box(15) .ne. 0 ) then
        if ( box(15) .ne. 0 ) then
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist1,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist1)
c
        nlist = nlist1
c
c       ... evaluate self interactions
c
         if ( nkids .eq. 0 ) then
         list(nlist+1) = ibox
         nlist = nlist + 1 
         endif
c    
c
c       ... evaluate interactions with list #3 
c
        itype=3
        call d3tgetl(ier,ibox,itype,list(nlist+1),nlist3,wlists)
        if (ifprint .ge. 2) call prinf('list3=*',list,nlist3)
        nlist = nlist + nlist3
c
c    
c       ... evaluate interactions with list #4
c
        itype=4
        call d3tgetl(ier,ibox,itype,list(nlist+1),nlist4,wlists)
        if (ifprint .ge. 2) call prinf('list4=*',list,nlist4)
        nlist = nlist + nlist4
c
ccc        write(*,*) ibox, nlist, nlist1, nlist3, nlist4
c
c
c       ... for all pairs in list #1, #3, and #4
c       evaluate the potentials and fields directly
c
            do ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            if( ibox .eq. jbox ) ifself = 1
            if( ibox .ne. jbox ) ifself = 0

            call elhfmm3dtriaiter_count(ifself,box,box1,
     $         ifsingle,ifdouble,ifptfrc,ifstrain,
     $         ifptfrctarg,ifstraintarg,lused)
            ntotal = ntotal + lused
            lstore = lstore + lused 
            lpairs = lpairs + 1
c
            iwpairs(1,next_pair) = ibox
            iwpairs(2,next_pair) = jbox
            iwpairs(3,next_pair) = next_chunk
            iwpairs(4,next_pair) = lused
            next_pair = next_pair + 1
            next_chunk = next_chunk + lused
c
ccc            write(*,*) ibox, jbox, lused
c
            enddo
        endif
c
        npairs = npairs + lpairs
ccc        write(*,*) ibox, lpairs
c
        if( lpairs .eq. 0 ) iwboxes(1,ibox)=0
        if( lpairs .ne. 0 ) iwboxes(1,ibox)=lptr
        iwboxes(2,ibox)=lpairs
        if( lstore .eq. 0 ) iwboxes(3,ibox)=0
        if( lstore .ne. 0 ) iwboxes(3,ibox)=iptr
        iwboxes(4,ibox)=lstore
        lptr=lptr+lpairs
        iptr=iptr+lstore
c
 6200   continue
c
c        call prinf('next_pair=*',next_pair,1)
c        call prinf('next_chunk=*',next_chunk,1)
c
c        call prinf('iwboxes=*',iwboxes,4*nboxes)
c        call prinf('iwpairs=*',iwpairs,4*npairs)
c
        call prinf('nboxes=*',nboxes,1)
        call prinf('npairs=*',npairs,1)
        call prinf('ntotal=*',ntotal,1)
        call prinf('ntotal(k)=*',ntotal/1000,1)
        call prinf('ntotal(M)=*',ntotal/1000000,1)
        call prinf('ntotal(G)=*',ntotal/1000000000,1)
c               
        lused = ntotal
c
c       ... storage
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
        if( icomp_type .eq. 4 ) goto 9000
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
C$OMP$PRIVATE(icomp,ipairs,istore,iavail,lused_istore) 
C$OMP$PRIVATE(ii,nlist1,nlist3,nlist4) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
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
        icomp = icomp_type
        if( icomp_type .eq. 0 ) icomp = 1
        if( icomp_type .eq. 5 ) icomp = 2
c       ... provided work array is too short, direct evaluation only
        if( lwsave .lt. lused ) icomp = 0
c
c       ... adjust list 3 and list 4 processing
c
ccc        if (nkids .eq. 0 ) then
        if ( 1 .eq. 1 ) then
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist1,wlists)
        if (ifprint .eq. 1) call prinf('list1=*',list,nlist1)
c
        nlist = nlist1
c
c
c       ... evaluate self interactions
c
        if (nkids .eq. 0 ) then
         list(nlist+1) = ibox
         nlist = nlist + 1 
        endif
c
c    
c       ... evaluate interactions with list #3 
c
        itype=3
        call d3tgetl(ier,ibox,itype,list(nlist+1),nlist3,wlists)
        if (ifprint .ge. 2) call prinf('list3=*',list,nlist3)
        nlist = nlist + nlist3
c
c    
c       ... evaluate interactions with list #4
c
        itype=4
        call d3tgetl(ier,ibox,itype,list(nlist+1),nlist4,wlists)
        if (ifprint .ge. 2) call prinf('list4=*',list,nlist4)
        nlist = nlist + nlist4
c
ccc        write(*,*) nlist, nlist3, nlist4
ccc        pause
c
c
c       ... for all pairs in list #1, #3, and #4
c       evaluate the potentials and fields directly
c    
            do 6203 ilist=1,nlist
c
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c
c          ipairs = iwboxes(1,ibox)
c          istore = iwpairs(3,ipairs+ilist-1)
c          iavail = iwpairs(4,ipairs+ilist-1)
c
c          ... search for allocated workspace, slow
c
          ipairs = iwboxes(1,jbox)
          do ii=1,iwboxes(2,jbox)
          if( iwpairs(1,ipairs+ii-1) .ne. jbox ) then
          endif
          if( iwpairs(2,ipairs+ii-1) .eq. ibox ) then
          istore = iwpairs(3,ipairs+ii-1)
          iavail = iwpairs(4,ipairs+ii-1)
          endif
          enddo
c
c          call prinf('ibox=*',ibox,1)
c          call prinf('jbox=*',jbox,1)
c          call prinf('ipairs=*',ipairs,1)
c          call prinf('ilist=*',ilist,1)
c          call prinf('nlist=*',nlist,1)
c          call prinf('ipairs ptr=*',ipairs+ilist-1,1)
c          call prinf('istore=*',istore,1)
c          call prinf('iavail=*',iavail,1)
c
               call elhfmm3dtriaiter_direct(ifself,box1,box,
     $            RLAM,RMU,TRIANGLE,TRINORM,SOURCE,
     $            ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $            ifptfrc,ptfrc,ifstrain,strain,
     $            target,ifptfrctarg,PTFRCtarg,
     $            ifstraintarg,STRAINtarg,
     $            icomp,wsave(istore),lused_istore)
c
c               if( iavail .lt. lused_istore ) then
c               write(*,*) nlist, nlist1, nlist3, nlist4
c               write(*,*) ibox, jbox, iavail, lused_istore
c               pause
c               endif
c
ccc        call prin2('wsave=*',wsave(istore),lused_istore)

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
 9000   continue
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
        subroutine elhfmm3dtriaiter_count(ifself,box,box1,
     $         ifsingle,ifdouble,ifptfrc,ifstrain,
     $         ifptfrctarg,ifstraintarg,lused)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        lused = 0
c
c       ... Self interactions are always included, assume normally
c       oriented dipole, i.e. 9 words per dipole strength vector 
c       (need 27 words for arbitrary oriented dipole). 
c
        nsingle = 3
        if(ifsingle .eq. 1 ) then
        if(ifptfrc .eq. 1)
     $     lused = lused + nsingle*3*box(15)*box1(15)
        if(ifstrain .eq. 1) 
     $     lused = lused + nsingle*9*box(15)*box1(15)
        if(ifptfrctarg .eq. 1)
     $     lused = lused + nsingle*3*box(15)*box1(17)
        if(ifstraintarg .eq. 1) 
     $     lused = lused + nsingle*9*box(15)*box1(17)
        endif

        ndouble = 3
        if(ifdouble .eq. 1 ) then
        if(ifptfrc .eq. 1)
     $     lused = lused + ndouble*3*box(15)*box1(15)
        if(ifstrain .eq. 1) 
     $     lused = lused + ndouble*9*box(15)*box1(15)
        if(ifptfrctarg .eq. 1)
     $     lused = lused + ndouble*3*box(15)*box1(17)
        if(ifstraintarg .eq. 1) 
     $     lused = lused + ndouble*9*box(15)*box1(17)
        endif
c
c       
c       allocate 2 real *8 storage elements for each complex *16 element
c       
ccc        lused = lused*2
c
        if( lused .lt. 0 ) lused = 0
c
ccc        call prinf('lused=*',lused,1)
c      
        return
        end
c
c
c
c
c
        subroutine elhfmm3dtriaiter_direct(ifself,box,box1,
     $     RLAM,RMU,TRIANGLE,TRINORM,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg,
     $     icomp,wsave,lused)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triangle(3,3,1),trinorm(3,1),source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension sourceim(3),targetim(3)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
        dimension ptfrc1(3),strain1(3,3)
c
        dimension sigma_sl0(3), sigma_dl0(3)
c
        real *8 wsave(1)
c
        lused = 0
c
c       ... sources
c
        ione=1
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
        do j=box1(14),box1(14)+box1(15)-1
        sourceim(1)=source(1,j)
        sourceim(2)=source(2,j)
        sourceim(3)=-source(3,j)

        do i=box(14),box(14)+box(15)-1
c
        if (ifsingle .eq. 1 ) then

        if( icomp .eq. 1 ) then

        call eluh3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl(1,i),trinorm(1,i),sourceim,
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)

        endif

        if( icomp .eq. 2 .or. icomp .eq. 3 ) then

        if( ifptfrc .eq. 1 ) call elhfmm3d_arrzero(ptfrc0,3)
        if( ifstrain .eq. 1 ) call elhfmm3d_arrzero(strain0,3*3)
        
        do k=1,3
c
        if( icomp .eq. 2 ) then
        sigma_sl0(1)=0
        sigma_sl0(2)=0
        sigma_sl0(3)=0
        sigma_sl0(k)=1

        call eluh3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl0,trinorm(1,i),sourceim,
     1     ifptfrc,ptfrc1,ifstrain,strain1,numfunev)

        if( ifptfrc .eq. 1 ) then
        call elhfmm3d_arrmove(ptfrc1,wsave(lused+1),3)
        lused=lused+3
        endif
        if( ifstrain .eq. 1 ) then
        call elhfmm3d_arrmove(strain1,wsave(lused+1),9)
        lused=lused+9
        endif
        endif
c
        if( icomp .eq. 3 ) then
        if( ifptfrc .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),ptfrc1,3)
        lused=lused+3
        endif
        if( ifstrain .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),strain1,9)
        lused=lused+9
        endif
        endif
c
        if( ifptfrc .eq. 1 ) 
     $     call elhfmm3d_arrcadd(ptfrc1,sigma_sl(k,i),ptfrc0,3)
        if( ifstrain .eq. 1 ) 
     $     call elhfmm3d_arrcadd(strain1,sigma_sl(k,i),strain0,3*3)

        enddo        
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

        if( icomp .eq. 1 ) then

        call elth3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),sourceim,
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)

        endif

        if( icomp .eq. 2 .or. icomp .eq. 3 ) then

        if( ifptfrc .eq. 1 ) call elhfmm3d_arrzero(ptfrc0,3)
        if( ifstrain .eq. 1 ) call elhfmm3d_arrzero(strain0,3*3)
        
        do k=1,3
c
        if( icomp .eq. 2 ) then
        sigma_dl0(1)=0
        sigma_dl0(2)=0
        sigma_dl0(3)=0
        sigma_dl0(k)=1

        call elth3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl0,sigma_dv(1,i),sourceim,
     1     ifptfrc,ptfrc1,ifstrain,strain1,numfunev)

        if( ifptfrc .eq. 1 ) then
        call elhfmm3d_arrmove(ptfrc1,wsave(lused+1),3)
        lused=lused+3
        endif
        if( ifstrain .eq. 1 ) then
        call elhfmm3d_arrmove(strain1,wsave(lused+1),9)
        lused=lused+9
        endif
        endif
c
        if( icomp .eq. 3 ) then
        if( ifptfrc .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),ptfrc1,3)
        lused=lused+3
        endif
        if( ifstrain .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),strain1,9)
        lused=lused+9
        endif
        endif
c
        if( ifptfrc .eq. 1 ) 
     $     call elhfmm3d_arrcadd(ptfrc1,sigma_dl(k,i),ptfrc0,3)
        if( ifstrain .eq. 1 ) 
     $     call elhfmm3d_arrcadd(strain1,sigma_dl(k,i),strain0,3*3)

        enddo        
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
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
        do j=box1(16),box1(16)+box1(17)-1
c
        targetim(1)=target(1,j)
        targetim(2)=target(2,j)
        targetim(3)=-target(3,j)

        do i=box(14),box(14)+box(15)-1
        if (ifsingle .eq. 1 ) then

        if( icomp .eq. 1 ) then

         call eluh3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl(1,i),trinorm(1,i),targetim,
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)

        endif

        if( icomp .eq. 2 .or. icomp .eq. 3 ) then

        if( ifptfrctarg .eq. 1 ) call elhfmm3d_arrzero(ptfrc0,3)
        if( ifstraintarg .eq. 1 ) call elhfmm3d_arrzero(strain0,3*3)
        
        do k=1,3
c
        if( icomp .eq. 2 ) then
        sigma_sl0(1)=0
        sigma_sl0(2)=0
        sigma_sl0(3)=0
        sigma_sl0(k)=1

        call eluh3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl0,trinorm(1,i),targetim,
     1     ifptfrctarg,ptfrc1,ifstraintarg,strain1,numfunev)

        if( ifptfrctarg .eq. 1 ) then
        call elhfmm3d_arrmove(ptfrc1,wsave(lused+1),3)
        lused=lused+3
        endif
        if( ifstraintarg .eq. 1 ) then
        call elhfmm3d_arrmove(strain1,wsave(lused+1),9)
        lused=lused+9
        endif
        endif
c
        if( icomp .eq. 3 ) then
        if( ifptfrctarg .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),ptfrc1,3)
        lused=lused+3
        endif
        if( ifstraintarg .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),strain1,9)
        lused=lused+9
        endif
        endif
c
        if( ifptfrctarg .eq. 1 ) 
     $     call elhfmm3d_arrcadd(ptfrc1,sigma_sl(k,i),ptfrc0,3)
        if( ifstraintarg .eq. 1 ) 
     $     call elhfmm3d_arrcadd(strain1,sigma_sl(k,i),strain0,3*3)

        enddo        
        endif

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

        if( icomp .eq. 1 ) then

        call elth3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),targetim,
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)

        endif

        if( icomp .eq. 2 .or. icomp .eq. 3 ) then

        if( ifptfrctarg .eq. 1 ) call elhfmm3d_arrzero(ptfrc0,3)
        if( ifstraintarg .eq. 1 ) call elhfmm3d_arrzero(strain0,3*3)
        
        do k=1,3
c
        if( icomp .eq. 2 ) then
        sigma_dl0(1)=0
        sigma_dl0(2)=0
        sigma_dl0(3)=0
        sigma_dl0(k)=1

        call elth3triaadap        
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl0,sigma_dv(1,i),targetim,
     1     ifptfrctarg,ptfrc1,ifstraintarg,strain1,numfunev)

        if( ifptfrctarg .eq. 1 ) then
        call elhfmm3d_arrmove(ptfrc1,wsave(lused+1),3)
        lused=lused+3
        endif
        if( ifstraintarg .eq. 1 ) then
        call elhfmm3d_arrmove(strain1,wsave(lused+1),9)
        lused=lused+9
        endif
        endif
c
        if( icomp .eq. 3 ) then
        if( ifptfrctarg .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),ptfrc1,3)
        lused=lused+3
        endif
        if( ifstraintarg .eq. 1 ) then
        call elhfmm3d_arrmove(wsave(lused+1),strain1,9)
        lused=lused+9
        endif
        endif
c
        if( ifptfrctarg .eq. 1 ) 
     $     call elhfmm3d_arrcadd(ptfrc1,sigma_dl(k,i),ptfrc0,3)
        if( ifstraintarg .eq. 1 ) 
     $     call elhfmm3d_arrcadd(strain1,sigma_dl(k,i),strain0,3*3)

        enddo        
        endif


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
        subroutine elhfmm3d_arrmove(x,y,n)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1)
c
        do 1200 i=1,n
        y(i)=x(i)
 1200   continue
c
        return
        end
c
c
c
c
c
        subroutine elhfmm3d_arrzero(y,n)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1)
c
        do 1200 i=1,n
        y(i)=0
 1200   continue
c
        return
        end
c
c
c
c
c
        subroutine elhfmm3d_arrcadd(x,c,y,n)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1)
c
        do 1200 i=1,n
        y(i)=y(i)+c*x(i)
 1200   continue
c
        return
        end
c
c
c
c
c
        subroutine eluh3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,ifptfrc,ptfrc,ifstrain,strain,numfunev)
C
C     adaptive integration on triangles for elastostatic single layer.
C
C     INPUT:
C
C     RLAM, RMU = Lame parameters
C     NTRI = number of triangles
C     TRIANGLES(3,3,ntri) = (x,y,z)-coordinate of 3 vertices 
C     SIGMA(3,ntri) = vector strength of nth charge
C     TRINORM(3,ntri) =  normal vectors
C     TARGET(3) = evaluation point
C
C     OUTPUT:
C
C     PTFRC(3) = computed displacement at target
C     STRAIN(3,3) = computed strain at target
C
        implicit real *8 (a-h,o-z)
        real *8 triangles(3,3,ntri),sigma(3,ntri)
        real *8 trinorm(3,ntri)
        real *8 target(3),source(3),ptfrc(3),strain(3,3)
        real *8 ptfrc0(3),strain0(3,3)
        real *8 vert1(3),vert2(3),vert3(3)
        dimension rints(100),par(100)
	real *8, allocatable :: w(:)
        integer ntri
c
        external fun3eluh_eval
c
	allocate( w(100000) )
c
        ptfrc(1) = 0.0d0
        ptfrc(2) = 0.0d0
        ptfrc(3) = 0.0d0
c
        do i=1,3
        do j=1,3
        strain(i,j) = 0.0d0
        enddo
        enddo
c       
c       ... note that the sum of weights is .5, 
c       we will need to include a factor of 2 later
c
        numfunev=0
c
        do 1200 itri = 1,ntri
c
        vert1(1) = triangles(1,1,itri)
        vert1(2) = triangles(2,1,itri)
        vert1(3) = triangles(3,1,itri)
        vert2(1) = triangles(1,2,itri)
        vert2(2) = triangles(2,2,itri)
        vert2(3) = triangles(3,2,itri)
        vert3(1) = triangles(1,3,itri)
        vert3(2) = triangles(2,3,itri)
        vert3(3) = triangles(3,3,itri)
        
        par(1)=rlam
        par(2)=rmu
        par(3)=sigma(1,itri)
        par(4)=sigma(2,itri)
        par(5)=sigma(3,itri)
        par(6)=trinorm(1,itri)
        par(7)=trinorm(2,itri)
        par(8)=trinorm(3,itri)
        par(9)=ifptfrc
        par(10)=ifstrain

        nq = 6
        nfuns = 12
c       
        eps=1e-6
ccc        eps=1e-3
        call tria3adam(ier,vert1,vert2,vert3,fun3eluh_eval,nfuns,
     1      target,par,nq,eps,rints,maxrec,numfunev0,w)
ccc        call prinf('ier=*',ier,1)
ccc        call prinf('numfunev=*',numfunev0,1)
        numfunev=numfunev+numfunev0

        k=0
        
        if( ifptfrc .eq. 1 ) then
        do i = 1,3
        k=k+1
        ptfrc0(i)=rints(k)
        enddo
        endif
        
        if( ifstrain .eq. 1 ) then
        do i = 1,3
        do j = 1,3
        k=k+1
        strain0(i,j)=rints(k)
        enddo
        enddo
        endif
c
c
        if( ifptfrc .eq. 1 ) then
        do i = 1,3
        ptfrc(i) = ptfrc(i) + ptfrc0(i) 
        enddo
        endif
        
        if( ifstrain .eq. 1 ) then
        do i = 1,3
        do j = 1,3
        strain(i,j) = strain(i,j) + strain0(i,j) 
        enddo
        enddo
        endif
c
 1200   continue
        return
        end
c
c
c
c
c
        subroutine elth3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,ifptfrc,ptfrc,ifstrain,strain,numfunev)
C
C     adaptive integration on triangles for elastostatic double layer.
C
C     INPUT:
C
C     RLAM, RMU = Lame parameters
C     NTRI = number of triangles
C     TRIANGLES(3,3,ntri) = (x,y,z)-coordinate of 3 vertices 
C     SIGMA(3,ntri) = vector strength of nth charge
C     TRINORM(3,ntri) =  normal vectors
C     TARGET(3) = evaluation point
C
C     OUTPUT:
C
C     PTFRC(3) = computed displacement at target
C     STRAIN(3,3) = computed strain at target
C
        implicit real *8 (a-h,o-z)
        real *8 triangles(3,3,ntri),sigma(3,ntri)
        real *8 trinorm(3,ntri)
        real *8 target(3),source(3),ptfrc(3),strain(3,3)
        real *8 ptfrc0(3),strain0(3,3)
        real *8 vert1(3),vert2(3),vert3(3)
        dimension rints(100),par(100)
	real *8, allocatable :: w(:)
        integer ntri
c
        external fun3elth_eval
c
	allocate( w(100000) )
c
        ptfrc(1) = 0.0d0
        ptfrc(2) = 0.0d0
        ptfrc(3) = 0.0d0
c
        do i=1,3
        do j=1,3
        strain(i,j) = 0.0d0
        enddo
        enddo
c       
c       ... note that the sum of weights is .5, 
c       we will need to include a factor of 2 later
c
        numfunev=0
c
        do 1200 itri = 1,ntri
c
        vert1(1) = triangles(1,1,itri)
        vert1(2) = triangles(2,1,itri)
        vert1(3) = triangles(3,1,itri)
        vert2(1) = triangles(1,2,itri)
        vert2(2) = triangles(2,2,itri)
        vert2(3) = triangles(3,2,itri)
        vert3(1) = triangles(1,3,itri)
        vert3(2) = triangles(2,3,itri)
        vert3(3) = triangles(3,3,itri)
        
        par(1)=rlam
        par(2)=rmu
        par(3)=sigma(1,itri)
        par(4)=sigma(2,itri)
        par(5)=sigma(3,itri)
        par(6)=trinorm(1,itri)
        par(7)=trinorm(2,itri)
        par(8)=trinorm(3,itri)
        par(9)=ifptfrc
        par(10)=ifstrain
        
        nq = 6
        nfuns = 12
c       
        eps=1e-6
ccc        eps=1e-3
        call tria3adam(ier,vert1,vert2,vert3,fun3elth_eval,nfuns,
     1      target,par,nq,eps,rints,maxrec,numfunev0,w)
ccc        call prinf('ier=*',ier,1)
ccc        call prinf('numfunev=*',numfunev0,1)
        numfunev=numfunev+numfunev0

        k=0

        if( ifptfrc .eq. 1 ) then 
        do i = 1,3
        k=k+1
        ptfrc0(i)=rints(k)
        enddo
        endif

        if( ifstrain .eq. 1 ) then
        do i = 1,3
        do j = 1,3
        k=k+1
        strain0(i,j)=rints(k)
        enddo
        enddo
        endif
c
c
        if( ifptfrc .eq. 1 ) then
        do i = 1,3
        ptfrc(i) = ptfrc(i) + ptfrc0(i) 
        enddo
        endif
        
        if( ifstrain .eq. 1 ) then
        do i = 1,3
        do j = 1,3
        strain(i,j) = strain(i,j) + strain0(i,j) 
        enddo
        enddo
        endif
c
 1200   continue
        return
        end
c
c
c
c
c
        subroutine fun3eluh_eval(x,y,z,target,par,f)
        implicit real *8 (a-h,o-z)
        dimension target(3),source(3),rvec(3)
        dimension ptfrc0(3),strain0(3,3),f(*)
        dimension sigma(3),trinorm(3),par(*)
c
        source(1)=x
        source(2)=y
        source(3)=z
c
        rlam=par(1)
        rmu=par(2)
        sigma(1)=par(3)
        sigma(2)=par(4)
        sigma(3)=par(5)
        trinorm(1)=par(6)
        trinorm(2)=par(7)
        trinorm(3)=par(8)
c
        ifptfrc = par(9)
        ifstrain = par(10)
c
c        call green3elu_eval
c     $     (rlam,rmu,source,sigma,target,
c     $     ptfrc0,ifstrain,strain0)
c
c        call green3eluh_eval
c     $     (rlam,rmu,source,sigma,target,
c     $     ptfrc0,ifstrain,strain0)
c
        call green3eluh_image_eval
     $     (rlam,rmu,source,sigma,target,
     $     ptfrc0,ifstrain,strain0)
c
c        call green3eluh_image_a_eval
c     $     (rlam,rmu,source,sigma,target,
c     $     ptfrc0,ifstrain,strain0)
c
c        call green3eluh_mindlin_eval
c     $     (rlam,rmu,source,sigma,target,
c     $     ptfrc0,ifstrain,strain0)
c
        k=0
        
        if( ifptfrc .eq. 1) then
        do i = 1,3
        k=k+1
        f(k) = ptfrc0(i)
        enddo
        endif

        if( ifstrain .eq. 1 ) then
        do i = 1,3
        do j = 1,3
        k=k+1
        f(k) = strain0(i,j)
        enddo
        enddo
        endif
c
        return
        end
c
c
c
        subroutine fun3elth_eval(x,y,z,target,par,f)
        implicit real *8 (a-h,o-z)
        dimension target(3),source(3),rvec(3)
        dimension ptfrc0(3),strain0(3,3),f(*)
        dimension sigma(3),trinorm(3),par(*)
c
        source(1)=x
        source(2)=y
        source(3)=z
c
        rlam=par(1)
        rmu=par(2)
        sigma(1)=par(3)
        sigma(2)=par(4)
        sigma(3)=par(5)
        trinorm(1)=par(6)
        trinorm(2)=par(7)
        trinorm(3)=par(8)
c
        ifptfrc = par(9)
        ifstrain = par(10)
c
c        call green3elt_eval
c     $     (rlam,rmu,source,sigma,trinorm,target,
c     $     ptfrc0,ifstrain,strain0)
c
c        call green3elth_eval
c     $     (rlam,rmu,source,sigma,trinorm,target,
c     $     ptfrc0,ifstrain,strain0)
c
        call green3elth_image_eval
     $     (rlam,rmu,source,sigma,trinorm,target,
     $     ptfrc0,ifstrain,strain0)
c
c        call green3elth_image_a_eval
c     $     (rlam,rmu,source,sigma,trinorm,target,
c     $     ptfrc0,ifstrain,strain0)
c
c        call green3elth_mindlin_eval
c     $     (rlam,rmu,source,sigma,trinorm,target,
c     $     ptfrc0,ifstrain,strain0)
c
        k=0
c
        if( ifptfrc .eq. 1 ) then
        do i = 1,3
        k=k+1
        f(k) = ptfrc0(i)
        enddo
        endif
 
        if( ifstrain .eq. 1 ) then
        do i = 1,3
        do j = 1,3
        k=k+1
        f(k) = strain0(i,j)
        enddo
        enddo
        endif
c
        return
        end
c
c
c
c
c
        subroutine elh3dtriadirecttarg(
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     NTARGET,target,
     $     ifptfrctarg,PTFRCtarg,ifstraintarg,STRAINtarg)
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
C$OMP$PRIVATE(i,j,ptfrc0,strain0,numfunev)
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
c
c       ... image contribution
c
        if (ifsingle .eq. 1 ) then
        call eluh3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,trinorm,source(1,j),
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if( ifstrain .eq. 1 ) then
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

        if (ifdouble .eq. 1 ) then
        call elth3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,source(1,j),
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if( ifstrain .eq. 1 ) then
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
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0,numfunev)
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
c
c       ... image contribution
c
        if (ifsingle .eq. 1 ) then
        call eluh3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,trinorm,target(1,j),
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)
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

        if (ifdouble .eq. 1 ) then
        call elth3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,target(1,j),
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)
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
