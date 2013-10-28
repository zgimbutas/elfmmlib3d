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
      subroutine l3dtriaformmp_mindlinb(nqtri,scale,rlame,
     1           triangle,trianorm,source,ifsingle,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,center,nterms,mpwork,chwork,quadvec,mpole)
C***********************************************************************
C
C     This subroutines forms the multipole expansion due to layer potential
C     densities on a collection of triangles,
C
C     INPUT:
C
C     nqtri       : order of desired quadrature rule
C     scale       : scaling parameter for local expansion
C     rlame       : Lame coefficients
C     triangle    : array of triangle coordinates
C     trianorm    : array of triangle normals
C     source      : array of triangle centroids
C     ifsingle    : flag for single layer source (1 => present)
C     charge      : single layer force vector
C     ifdouble    : flag for double layer source (1 => present)
C     dipstr      : double layer force vector
C     dipvec      : double layer normal vector
C     ns          : number of sources
C     center      : expansion center
C     nterms      : order of spherical harmonic expansion
C     mpwork      : work array
C     chwork      : work array
C     quadvec     : work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier         : error code for memory allocation inside routine
C     local       : coefficients of s.h. expansion
C---------------------------------------------------------------------
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
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpwork(0:nterms,-nterms:nterms)
      complex *16 chwork(ns)
      real *8 quadvec(6,ns)
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
      call formmp_mindlinb(scale,rlame,source_interp,ifsingle,
     1           charge_interp,ifdouble,dipstr_interp,dipvec_interp,
     1           nsource,center,nterms,mpwork,chwork,quadvec,mpole)
c
      return
      end
c
c
c
c
      subroutine l3dtriaformmp_mindlinc(nqtri,scale,rlame,
     1           triangle,trianorm,source,ifsingle,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,center,nterms,mpwork,chwork,quadvec,
     1           dipvec2,mpole,mpole2)
C***********************************************************************
C
C     This subroutines forms the Mindlin C mpole expansions due to 
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
C     nterms      : order of spherical harmonic expansion
C     mpwork      : work array
C     chwork      : work array
C     quadvec     : work array
C     dipvec2     : work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier         : error code for memory allocation inside routine
C     mpole       : coefficients of first Mindlin C s.h. expansion
C     mpole2      : coefficients of second Mindlin C  s.h. expansion
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
      real *8 quadvec(6,ns)
      real *8 dipvec2(3,ns)
      complex *16 mpwork(0:nterms,-nterms:nterms)
      complex *16 chwork(ns)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
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
      call formmp_mindlinc(scale,rlame,source_interp,ifsingle,
     1           charge_interp,ifdouble,dipstr_interp,dipvec_interp,
     1           nsource,center,nterms,mpwork,chwork,quadvec,dipvec2,
     1           mpole,mpole2)
c
      return
      end





      subroutine l3dtriaformmp_mindlinb_trunc(nqtri,scale,rlame,
     1           triangle,trianorm,source,ifsingle,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,center,nterms,
     $           mpole,wlege,nlege)
C***********************************************************************
C
C     This subroutines forms the multipole expansion due to layer potential
C     densities on a collection of triangles,
C
C     INPUT:
C
C     nqtri       : order of desired quadrature rule
C     scale       : scaling parameter for local expansion
C     rlame       : Lame coefficients
C     triangle    : array of triangle coordinates
C     trianorm    : array of triangle normals
C     source      : array of triangle centroids
C     ifsingle    : flag for single layer source (1 => present)
C     charge      : single layer force vector
C     ifdouble    : flag for double layer source (1 => present)
C     dipstr      : double layer force vector
C     dipvec      : double layer normal vector
C     ns          : number of sources
C     center      : expansion center
C     nterms      : order of spherical harmonic expansion
C     mpwork      : work array
C     chwork      : work array
C     quadvec     : work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier         : error code for memory allocation inside routine
C     local       : coefficients of s.h. expansion
C---------------------------------------------------------------------
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
      complex *16 mpole(0:nterms,-nterms:nterms)
ccc      complex *16 mpwork(0:nterms,-nterms:nterms)
ccc      complex *16 chwork(ns)
ccc      real *8 quadvec(6,ns)
      complex *16, allocatable :: mpwork(:,:)
      complex *16, allocatable :: chwork(:)
      real *8, allocatable :: quadvec(:,:)
      integer nlege
      real *8 wlege(*)
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
        allocate(mpwork(0:nterms,-nterms:nterms))        
        allocate(chwork(nsource))
        allocate(quadvec(6,nsource))
c
      call formmp_mindlinb_trunc(scale,rlame,source_interp,ifsingle,
     1           charge_interp,ifdouble,dipstr_interp,dipvec_interp,
     1           nsource,center,nterms,mpwork,chwork,quadvec,mpole,
     $           wlege,nlege)
c
      return
      end
c
c
c
c
      subroutine l3dtriaformmp_mindlinc_trunc(nqtri,scale,rlame,
     1           triangle,trianorm,source,ifsingle,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,center,nterms,
     1           mpole,mpole2,wlege,nlege)
C***********************************************************************
C
C     This subroutines forms the Mindlin C mpole expansions due to 
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
C     nterms      : order of spherical harmonic expansion
C     mpwork      : work array
C     chwork      : work array
C     quadvec     : work array
C     dipvec2     : work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier         : error code for memory allocation inside routine
C     mpole       : coefficients of first Mindlin C s.h. expansion
C     mpole2      : coefficients of second Mindlin C  s.h. expansion
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
ccc      complex *16 mpwork(0:nterms,-nterms:nterms)
ccc      complex *16 chwork(ns)
ccc      real *8 quadvec(6,ns)
ccc      real *8 dipvec2(3,ns)
      complex *16, allocatable :: mpwork(:,:)
      complex *16, allocatable :: chwork(:)
      real *8, allocatable :: quadvec(:,:)
      real *8, allocatable :: dipvec2(:,:)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
      integer nlege
      real *8 wlege(*)
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
        allocate(mpwork(0:nterms,-nterms:nterms))        
        allocate(chwork(nsource))
        allocate(quadvec(6,nsource))
        allocate(dipvec2(3,nsource))
c
      call formmp_mindlinc_trunc(scale,rlame,source_interp,ifsingle,
     1           charge_interp,ifdouble,dipstr_interp,dipvec_interp,
     1           nsource,center,nterms,mpwork,chwork,quadvec,dipvec2,
     1           mpole,mpole2,wlege,nlege)
c
      return
      end
