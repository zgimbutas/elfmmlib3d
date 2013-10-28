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
c        This file contains the routines for elastostatic layer
c        potentials in half space in R^3. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      elh3dtriadirecttarg - evaluates the elastostatic potential ON OR
c         OFF SURFACE due to a collection of flat triangles with
c         constant single and/or double layer densities using the direct
c         O(N^2) algorithm. Half space.
c
c      eluh3triaadap - evaluates the elastostatic potential for the
c         single layer Mindlin image (A+B+C) contribution at an
c         abritrary target due to a collection of flat triangles with
c         piecewise constant single layer density VIA DIRECT ADAPTIVE
c         QUADRATURE.
c
c      elth3triaadap - evaluates the elastostatic potential for the
c         double layer Mindlin image (A+B+C) contribution at an
c         abritrary target due to a collection of flat triangles with
c         piecewise constant single layer density VIA DIRECT ADAPTIVE
c         QUADRATURE.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c       O(N^2) algorithm. Half space.
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
        subroutine eluh3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,ifptfrc,ptfrc,ifstrain,strain,numfunev)
C
C     Adaptive integration on triangles for elastostatic single layer.
C     Half space Mindlin image (A+B+C).
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

        nq = 20
        nfuns = 12
c       
        eps=1e-12
ccc        eps=1e-6
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
C     Adaptive integration on triangles for elastostatic double layer.
C     Half space Mindlin image (A+B+C).
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
        
        nq = 20
        nfuns = 12
c       
        eps=1e-12
ccc        eps=1e-6
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
