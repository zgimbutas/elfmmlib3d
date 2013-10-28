cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    Subroutines for computing local expansions due to triangles.
c
c    l3dformtaintkerbf90
c    l3dformtaintkercf90
c
c
c
c
c
c
      subroutine l3dtriaformtaintkerbf90(ier,nqtri,scale,
     1           triangle,trianorm,source,ifsingle,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,rlame,center,radius,nterms,local)
C***********************************************************************
C
C     This subroutines forms the local expansion due to layer potential
C     densities on a collection of triangles,
C
C     INPUT:
C
C     nqtri       : order of desired quadrature rule
C     scale       : scaling parameter for local expansion
C     triangle    : array of triangle coordinates
C     trianorm    : array of triangle normals
C     source      : array of triangle centroids
C     ifsingle    : flag for single layer source (1 => present)
C     charge      : single layer force vector
C     ifdouble    : flag for double layer source (1 => present)
C     dipstr      : double layer force vector
C     dipvec      : double layer normal vector
C     ns          : number of sources
C     rlame       : Lame coefficients
C     center      : expansion center
C     radius      : sphere radius
C     nterms      : order of spherical harmonic expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier         : error code for memory allocation inside routine
C     local       : coefficients of s.h. expansion
C---------------------------------------------------------------------

C
      implicit real *8 (a-h,o-z)
      integer ier,iflag,lw,lused,ns,nquad,nquadm,nterms
      real *8 triangle(3,3,ns)
      real *8 trianorm(3,ns)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      real *8 source(3,ns),center(3)
      real *8 rlame(2)
      real *8 vert1(3),vert2(3),vert3(3)
      real *8 vertout(3)
      real *8 vert01(3),vert02(3),vert03(3)
      real *8 rnodes(2,1000)
      real *8 weights(1000)
      complex *16 local(0:nterms,-nterms:nterms)
c
      real *8, allocatable :: source_interp(:,:,:)
      real *8, allocatable :: charge_interp(:,:,:)
      real *8, allocatable :: dipstr_interp(:,:,:)
      real *8, allocatable :: dipvec_interp(:,:,:)
c
c     get smooth quadrature nodes/weights for each triangle and 
c     call point interaction
c
      vert01(1)=0
      vert01(2)=0
      vert02(1)=1
      vert02(2)=0
      vert03(1)=0
      vert03(2)=1
      call triasymq(nqtri,vert01,vert02,vert03,rnodes,weights,nnodes)
cc        call prinf('after triasymq, nqtri=*',nqtri,1)
cc        call prinf('after triasymq, nnodes=*',nnodes,1)
c        call prin2('after triasymq, rnodes=*',rnodes,2*nnodes)
c        call prin2('after triasymq, weights=*',weights,nnodes)
c
      ntri = ns
ccc      call prinf( ' ntri is *',ntri,1)
ccc      call prinf( ' nqtri is *',nqtri,1)
c
      allocate( source_interp(3,nnodes,ntri) )
c
      if( ifsingle .eq. 1 ) then
         allocate( charge_interp(3,nnodes,ntri) )
      else
         allocate( charge_interp(3,nnodes,1) )
      endif
        
      if( ifdouble .eq. 1 ) then
         allocate( dipstr_interp(3,nnodes,ntri) )
         allocate( dipvec_interp(3,nnodes,ntri) )
      else
         allocate( dipstr_interp(3,nnodes,1) )
         allocate( dipvec_interp(3,nnodes,1) )
      endif
c
      nsource = ntri*nnodes
      do j=1,ntri
         call triangle_area(triangle(1,1,j),ds)
         do i=1,nnodes
            u=rnodes(1,i)
            v=rnodes(2,i)
            source_interp(1,i,j)=triangle(1,1,j)+u*(triangle(1,2,j)
     $        -triangle(1,1,j))+v*(triangle(1,3,j)-triangle(1,1,j))
            source_interp(2,i,j)=triangle(2,1,j)+u*(triangle(2,2,j)
     $        -triangle(2,1,j))+v*(triangle(2,3,j)-triangle(2,1,j))
            source_interp(3,i,j)=triangle(3,1,j)+u*(triangle(3,2,j)
     $        -triangle(3,1,j))+v*(triangle(3,3,j)-triangle(3,1,j))
        
            if( ifsingle .eq. 1 ) then
               charge_interp(1,i,j)=charge(1,j)*weights(i)*ds*2
               charge_interp(2,i,j)=charge(2,j)*weights(i)*ds*2
               charge_interp(3,i,j)=charge(3,j)*weights(i)*ds*2
            endif
            if( ifdouble .eq. 1 ) then
               dipstr_interp(1,i,j)=dipstr(1,j)*weights(i)*ds*2
               dipstr_interp(2,i,j)=dipstr(2,j)*weights(i)*ds*2
               dipstr_interp(3,i,j)=dipstr(3,j)*weights(i)*ds*2
               dipvec_interp(1,i,j)=dipvec(1,j)
               dipvec_interp(2,i,j)=dipvec(2,j)
               dipvec_interp(3,i,j)=dipvec(3,j)
            endif
         enddo
      enddo
c
      call l3dformtaintkerbf90(ier,scale,source_interp,ifsingle,
     1           charge_interp,ifdouble,dipstr_interp,dipvec_interp,
     1           nsource,rlame,center,radius,nterms,local)
c
      return
      end
c
c
c
c
      subroutine l3dtriaformtaintkercf90(ier,nqtri,scale,
     1           triangle,trianorm,source,ifsingle,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,rlame,center,radius,nterms,local,local2)
C***********************************************************************
C
C     This subroutines forms the Mindlin C local expansions due to 
C     layer potential densities on a collection of triangles,
C
C     INPUT:
C
C     nqtri       : order of desired quadrature rule
C     scale       : scaling parameter for local expansion
C     triangle    : array of triangle coordinates
C     trianorm    : array of triangle normals
C     source      : array of triangle centroids
C     ifsingle    : flag for single layer source (1 => present)
C     charge      : single layer force vector
C     ifdouble    : flag for double layer source (1 => present)
C     dipstr      : double layer force vector
C     dipvec      : double layer normal vector
C     ns          : number of sources
C     rlame       : Lame coefficients
C     center      : expansion center
C     radius      : sphere radius
C     nterms      : order of spherical harmonic expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier         : error code for memory allocation inside routine
C     local       : coefficients of first Mindlin C s.h. expansion
C     local2       : coefficients of second Mindlin C  s.h. expansion
C---------------------------------------------------------------------

C
      implicit real *8 (a-h,o-z)
      integer ier,iflag,lw,lused,ns,nquad,nquadm,nterms
      real *8 triangle(3,3,ns)
      real *8 trianorm(3,ns)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      real *8 source(3,ns),center(3)
      real *8 rlame(2)
      real *8 vert1(3),vert2(3),vert3(3)
      real *8 vertout(3)
      real *8 vert01(3),vert02(3),vert03(3)
      real *8 rnodes(2,1000)
      real *8 weights(1000)
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 local2(0:nterms,-nterms:nterms)
c
      real *8, allocatable :: source_interp(:,:,:)
      real *8, allocatable :: charge_interp(:,:,:)
      real *8, allocatable :: dipstr_interp(:,:,:)
      real *8, allocatable :: dipvec_interp(:,:,:)
c
c     get smooth quadrature nodes/weights for each triangle and 
c     call point interaction
c
      vert01(1)=0
      vert01(2)=0
      vert02(1)=1
      vert02(2)=0
      vert03(1)=0
      vert03(2)=1
      call triasymq(nqtri,vert01,vert02,vert03,rnodes,weights,nnodes)
cc        call prinf('after triasymq, nqtri=*',nqtri,1)
cc        call prinf('after triasymq, nnodes=*',nnodes,1)
c        call prin2('after triasymq, rnodes=*',rnodes,2*nnodes)
c        call prin2('after triasymq, weights=*',weights,nnodes)
c
      ntri = ns
ccc      call prinf( ' ntri is *',ntri,1)
ccc      call prinf( ' nqtri is *',nqtri,1)
c
      allocate( source_interp(3,nnodes,ntri) )
c
      if( ifsingle .eq. 1 ) then
         allocate( charge_interp(3,nnodes,ntri) )
      else
         allocate( charge_interp(3,nnodes,1) )
      endif
        
      if( ifdouble .eq. 1 ) then
         allocate( dipstr_interp(3,nnodes,ntri) )
         allocate( dipvec_interp(3,nnodes,ntri) )
      else
         allocate( dipstr_interp(3,nnodes,1) )
         allocate( dipvec_interp(3,nnodes,1) )
      endif
c
      nsource = ntri*nnodes
      do j=1,ntri
         call triangle_area(triangle(1,1,j),ds)
         do i=1,nnodes
            u=rnodes(1,i)
            v=rnodes(2,i)
            source_interp(1,i,j)=triangle(1,1,j)+u*(triangle(1,2,j)
     $        -triangle(1,1,j))+v*(triangle(1,3,j)-triangle(1,1,j))
            source_interp(2,i,j)=triangle(2,1,j)+u*(triangle(2,2,j)
     $        -triangle(2,1,j))+v*(triangle(2,3,j)-triangle(2,1,j))
            source_interp(3,i,j)=triangle(3,1,j)+u*(triangle(3,2,j)
     $        -triangle(3,1,j))+v*(triangle(3,3,j)-triangle(3,1,j))
        
            if( ifsingle .eq. 1 ) then
               charge_interp(1,i,j)=charge(1,j)*weights(i)*ds*2
               charge_interp(2,i,j)=charge(2,j)*weights(i)*ds*2
               charge_interp(3,i,j)=charge(3,j)*weights(i)*ds*2
            endif
            if( ifdouble .eq. 1 ) then
               dipstr_interp(1,i,j)=dipstr(1,j)*weights(i)*ds*2
               dipstr_interp(2,i,j)=dipstr(2,j)*weights(i)*ds*2
               dipstr_interp(3,i,j)=dipstr(3,j)*weights(i)*ds*2
               dipvec_interp(1,i,j)=dipvec(1,j)
               dipvec_interp(2,i,j)=dipvec(2,j)
               dipvec_interp(3,i,j)=dipvec(3,j)
            endif
         enddo
      enddo
c
      call l3dformtaintkercf90(ier,scale,source_interp,ifsingle,
     1           charge_interp,ifdouble,dipstr_interp,dipvec_interp,
     1           nsource,rlame,center,radius,nterms,local,local2)
c
      return
      end

