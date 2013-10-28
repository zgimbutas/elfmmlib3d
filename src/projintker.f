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
c    Subroutines for projection onto local harmonic expansion
c
c
C***********************************************************************
      subroutine l3drescale(nterms,local,radius,scale)
C***********************************************************************
C
C     This subroutine rescales a spherical harmonic expansion
C     on the surface by r^n and scale^n, converting 
C     a surface function to the corresponding (scaled) harmonic 
C     expansion in the interior. Also rescale expansion coefficients 
C     by sqrt(2n+1) to compensate for the fact that l3dtaeval scales
c     by its inverse...
c
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           local = coefficients of s.h. expansion
C           radius = sphere radius
C           scale  = scale parameter for local expansion.
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = rescaled by sqrt(2n+1)/(scale^n r^n)
C---------------------------------------------------------------------
      implicit none
      integer nterms,l,m
      real *8 zmul,zscal,radius,scale
      complex *16 local(0:nterms,-nterms:nterms)
C
      zmul = 1.0d0
      zscal = 1.0d0/(scale*radius)
      do l=0,nterms
         do m=-l,l
	    local(l,m) = local(l,m)*zmul
	    local(l,m) = local(l,m)*sqrt(2*l+1.0d0)
         enddo
	 zmul = zmul*zscal
      enddo
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine l3devalsphereintker(source,charge,ns,intker,iflag,
     1           phival,center,radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine tabulates a harmonic function defined by 
C     a user-provided subroutine and a collection of sources/charges
C     at quadrature nodes on a sphere of radius RADIUS,
C     centered at CENTER.
C
C     The interaction kernel must take the form:
C
C     subroutine intker(iflag,source,charge,targ,pot,fld)
C
C     iflag determines what scalar function is used in projection.
C     If (iflag.eq.0) then intker returns the kernel in pot
C     If (iflag.eq.1) then intker returns the x-deriv of kernel in pot
C     If (iflag.eq.2) then intker returns the y-deriv of kernel in pot
C     If (iflag.eq.3) then intker returns the z-deriv of kernel in pot
C
C     In each case, intker returns the gradient of the scalar function
C     defined above in fld (NOT USED IN THIS SUBROUTINE).
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source   : source locations
C     charge   : source strength
C     ns       : number of sources
C     intker   : external subroutine
C     iflag    : flag determines whether the potential
C                defined by intker or one of its derivatives
C		 is used in defining scalar function to project
C     center   : sphere center
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     nquad    : number of quadrature nodes in theta direction
C     nquadm   : number of quadrature nodes in phi direction
C
C                total number of nodes on target sphere is nquad*nquadm
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit none
      integer nquad,nquadm,jj,kk,ns,iflag,nterms,i
      real *8 source(3,ns),targ(3), center(3)
      real *8 radius,xnodes(nquad)
      real *8 cosphi,sinphi,ctheta,stheta,phi,pi
      complex *16 phival(nquad,nquadm),charge(ns)
      complex *16 pot,fld(3)
      external intker
C
      pi = 4.0d0*datan(1.0d0)
      do jj=1,nquad
      do kk=1,nquadm
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
	 phi = 2*pi*kk/nquadm
	 cosphi = dcos(phi)
	 sinphi = dsin(phi)
	 targ(1) = center(1) + radius*stheta*cosphi
	 targ(2) = center(2) + radius*stheta*sinphi
	 targ(3) = center(3) + radius*ctheta
         phival(jj,kk) = 0.0d0
         pot = 0.0d0
         do i = 1,ns
            call intker(iflag,source(1,i),charge(i),targ,pot,fld)
            phival(jj,kk) = phival(jj,kk) + pot
         enddo
      enddo
      enddo
      return
      end
C
C
C
C
C
      subroutine l3dformtaintker(ier,scale,source,charge,ns,intker,
     1           iflag,center,radius,nterms,local,w,lw,lused)
C***********************************************************************
C
      implicit none
      integer ier,iflag,lw,lused,lused2,ns,nquad,nquadm,nterms
      integer iphival,lphival,ixnodes,iwts,iynm,lynm,iused,ifinit
      real *8 source(3,ns),center(3)
      real *8 scale,radius
      real *8 w(lw)
      complex *16 charge(ns)
      complex *16 local(0:nterms,-nterms:nterms)
      external intker
c
c
c
ccc      call prinf(' in formtaintker lw is *',lw,1)
      nquad = 2*nterms
      nquadm = nquad
      iphival = 1
      lphival = 2*nquad*nquadm
      ixnodes = iphival + lphival
      iwts = ixnodes + nquad
      iynm = iwts + nquad
      lynm = (nterms+1)**2
      iused = iynm + lynm
      if (iused.gt.lw) then
         ier = 1
         return
      endif
C
C     evaluate on sphere, project onto local expansion and rescale.
C
      ifinit = 1
      call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
      call l3devalsphereintker(source,charge,ns,intker,iflag,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call h3dprojloc(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local,w(iynm),w(iused),lw-iused,lused2)
      lused = iused+lused2
c
      call l3drescale(nterms,local,radius,scale)
      return
      end
C
C
C
C
      subroutine l3dformtaintkerf90(ier,scale,source,charge,ns,intker,
     1           iflag,center,radius,nterms,local)
C***********************************************************************
C
      implicit none
      integer ier,iflag,lw,lused,ns,nquad,nquadm,nterms
      integer ifinit,iphival,lphival,ixnodes,iwts,iynm,lynm
      integer iprojloc,iused,lused2
      real *8 source(3,ns),center(3)
      real *8 scale,radius
      real *8, allocatable :: w(:)
      complex *16 charge(ns)
      complex *16 local(0:nterms,-nterms:nterms)
      external intker
c
c
      nquad = 2*nterms
      nquadm = nquad
      iphival = 1
      lphival = 2*nquad*nquadm
      ixnodes = iphival + lphival
      iwts = ixnodes + nquad
      iynm = iwts + nquad
      lynm = (nterms+1)**2
      iprojloc = iynm + lynm
      iused = iprojloc + 2*nquad*(2*nterms+1)+10
      allocate(w(iused),stat=ier)
      if (ier.ne.0) then
         ier = 1
         return
      endif
C
C     evaluate on sphere, project onto local expansion and rescale.
C
      ifinit = 1
      call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
      call l3devalsphereintker(source,charge,ns,intker,iflag,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call h3dprojloc(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local,w(iynm),w(iprojloc),iused-iprojloc,lused2)
      lused = iused+lused2
c
      call l3drescale(nterms,local,radius,scale)
      return
      end
C
C
C
      subroutine l3dformtaintkervf90(ier,scale,source,ifcharge,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,intker,icomp,rlame,center,radius,nterms,local)
C***********************************************************************
C
      implicit none
      integer ier,iflag,lw,lused,ns,nquad,nquadm,nterms
      integer ifinit,iphival,lphival,ixnodes,iwts,iynm,lynm
      integer ifcharge,ifdouble,icomp,lused2,iprojloc,iused
      real *8 source(3,ns),center(3)
      real *8 scale,radius,rlame(2)
      real *8, allocatable :: w(:)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      complex *16 local(0:nterms,-nterms:nterms)
      external intker
c
c
      nquad = 2*nterms
      nquadm = nquad
      iphival = 1
      lphival = 2*nquad*nquadm
      ixnodes = iphival + lphival
      iwts = ixnodes + nquad
      iynm = iwts + nquad
      lynm = (nterms+1)**2
      iprojloc = iynm + lynm
      iused = iprojloc + 2*nquad*(2*nterms+1)+10
      allocate(w(iused),stat=ier)
      if (ier.ne.0) then
         ier = 1
         return
      endif
C
C     evaluate on sphere, project onto local expansion and rescale.
C
      ifinit = 1
      call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
      call l3devalsphereintkerv(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,intker,icomp,rlame,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call h3dprojloc(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local,w(iynm),w(iprojloc),iused-iprojloc,lused2)
      lused = iused+lused2
c
      call l3drescale(nterms,local,radius,scale)
      return
      end
c
c
c
C***********************************************************************
      subroutine l3devalsphereintkerv(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,intker,icomp,rlame,
     1     phival,center,radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine tabulates a harmonic function defined by 
C     a user-provided subroutine and a collection of sources/charges
C     at quadrature nodes on a sphere of radius RADIUS,
C     centered at CENTER.
C
C     The interaction kernel must take the form:
C
C     subroutine intker(icomp,rlame,source,ifcharge,charge,ifdouble,
C        dipstr,dipvec,targ,pot,iffld,fld)
C
C     icomp determines the component of interest
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source   : source locations
C     ifcharge : flag =1 if charge to be used
C     charge(3): source strength
C     ifdouble : flag =1 if dipole to be used
C     dipstr   : dipole strength
C     dipvec   : dipole vector
C     ns       : number of sources
C     intker   : external subroutine
C     icomp    : component of displacement
C     rlame    : Lame coefficients
C     center   : sphere center
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     nquad    : number of quadrature nodes in theta direction
C     nquadm   : number of quadrature nodes in phi direction
C
C                total number of nodes on target sphere is nquad*nquadm
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit none
      integer ns,ifcharge,ifdouble,icomp,nquad,nquadm,jj,kk
      integer nterms,i,iffld
      real *8 source(3,ns),targ(3), center(3), radius
      real *8 xnodes(nquad)
      real *8 rlame(2)
      real *8 cosphi,sinphi,ctheta,stheta,phi,pi
      complex *16 phival(nquad,nquadm)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      complex *16 pot,fld(3)
      external intker
C
ccc      call prinf('iflag is *',iflag,1)
      pi = 4.0d0*datan(1.0d0)
      do jj=1,nquad
      do kk=1,nquadm
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
	 phi = 2*pi*kk/nquadm
	 cosphi = dcos(phi)
	 sinphi = dsin(phi)
	 targ(1) = center(1) + radius*stheta*cosphi
	 targ(2) = center(2) + radius*stheta*sinphi
	 targ(3) = center(3) + radius*ctheta
         phival(jj,kk) = 0.0d0
         pot = 0.0d0
         iffld = 0
         do i = 1,ns
            call intker(icomp,rlame,source(1,i),ifcharge,charge(1,i),
     1           ifdouble,dipstr(1,i),dipvec(1,i),targ,pot,iffld,fld)
            phival(jj,kk) = phival(jj,kk) + pot
         enddo
      enddo
      enddo
      return
      end
C
