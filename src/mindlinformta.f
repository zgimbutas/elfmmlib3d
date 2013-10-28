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
c    l3drescale2
c    l3dformtaintkerbf90
c    l3devalsphereintkerb
c    intkerb
c    l3dformtaintkercf90
c    l3devalsphereintkerc1
c    intkerc1
c    l3devalsphereintkerc2
c    intkerc2
c
C***********************************************************************
      subroutine l3drescale2(nterms,local,radius,scfac)
C***********************************************************************
C
C     This subroutine rescales a spherical harmonic expansion
C     on the surface f radius radius by radius^n and scfac^n, 
C     converting a surface function to the corresponding (scaled) 
C     harmonic expansion in the interior. 
C     Also rescale expansion coefficients by sqrt(2n+1) to compensate 
C     for the fact that l3dtaeval scales by its inverse...
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nterms = order of spherical harmonic expansion
C     local = coefficients of s.h. expansion
C     radius = sphere radius
C     scfac  = scale parameter for local expansion.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = rescaled by sqrt(2n+1)/(scfac^n r^n)
C---------------------------------------------------------------------
      implicit none
      integer nterms,l,m
      real *8 radius
      real *8 scfac
      real *8 rmul,rscal
      complex *16 local(0:nterms,-nterms:nterms)
C
      rmul = 1.0d0
      rscal = 1.0d0/(scfac*radius)
      do l=0,nterms
         do m=-l,l
	    local(l,m) = local(l,m)*rmul
	    local(l,m) = local(l,m)*sqrt(2*l+1.0d0)
         enddo
	 rmul = rmul*rscal
      enddo
      return
      end
c
c
C***********************************************************************
      subroutine l3dformtaintkerbf90(ier,scfac,source,ifcharge,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,rlame,center,radius,nterms,local)
C***********************************************************************
C
C     This subrotine forms the local expansion due to a collection 
C     of Mindlin B sources.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     scfac     = scale parameter for local expansion.
C     source    = source locations
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ier    = error return code for memory allocation inside routine
C     local = coefficients of s.h. expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ier,lw,lused,ns,nquad,nquadm,nterms
      integer ifcharge,ifdouble
      integer ifinit
      integer iphival,lphival,ixnodes,iwts,iynm,lynm,iprojloc,iused
      real *8 source(3,ns),center(3)
      real *8 scfac,radius,rlame(2)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      complex *16 local(0:nterms,-nterms:nterms)
      real *8, allocatable :: w(:)
c
c     carve up workspace and compute allocation needed
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
      call l3devalsphereintkerb(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,rlame,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call projloc3d(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local,w(iprojloc),w(iynm))
c
      call l3drescale2(nterms,local,radius,scfac)
      return
      end
c
c
c
C***********************************************************************
      subroutine l3devalsphereintkerb(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,rlame,
     1     phival,center,radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine tabulates a harmonic function defined by 
C     a user-provided subroutine and a collection of sources/charges
C     at quadrature nodes on a sphere of radius RADIUS,
C     centered at CENTER.
C
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source    = source locations
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C     nquad     = number of quadrature nodes in theta direction
C     nquadm    = number of quadrature nodes in phi direction
C        [total number of nodes on target sphere is nquad*nquadm]
C     xnodes    = Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit none
      integer nterms,ns,ifcharge,ifdouble
      integer nquad,nquadm,jj,kk,i
      real *8 source(3,ns),targ(3), center(3)
      real *8 xnodes(nquad)
      real *8 rlame(2),radius
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      real *8 cosphi,ctheta,phi,pi,sinphi,stheta
      complex *16 phival(nquad,nquadm)
      complex *16 pot,fld(3)
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
            call intkerb(rlame,source(1,i),ifcharge,charge(1,i),
     1           ifdouble,dipstr(1,i),dipvec(1,i),targ,pot)
            phival(jj,kk) = phival(jj,kk) + pot
         enddo
      enddo
      enddo
      return
      end
C
C***********************************************************************
      subroutine intkerb(rlame,source,ifcharge,charge,
     1           ifdouble,dipstr,dipvec,targ,pot)
C***********************************************************************
C
C     This subroutine computes the pre-harmonic Mindlin B kernel
C     at a target point.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     rlame     = Lame coefficients
C     source    = source location
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     targ      = target location
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     pot       =  value of Mindlin B pre-harmonic potential 
C
C---------------------------------------------------------------------
      implicit none
      integer ifcharge,ifdouble
      integer nquad,nquadm,jj,kk
      real *8 source(3),targ(3)
      real *8 rlame(2)
      real *8 charge(3)
      real *8 dipstr(3)
      real *8 dipvec(3)
      real *8 alpha,scfac,r1,r2,r3,rr,d1,d2,pot2
      real *8 rnudotf,u11,u12,u22,u33
      complex *16 pot
C
      pot = 0.0d0
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      scfac  = (1.0d0-alpha)/alpha
      r1 = targ(1)-source(1)
      r2 = targ(2)-source(2)
      r3 = targ(3)-source(3)
      rr = sqrt(r1*r1+r2*r2+r3*r3)
      d1 = 1.0d0/(rr-r3)
      if (ifcharge.eq.1) then
         pot = pot + charge(1)*r1*d1
         pot = pot + charge(2)*r2*d1
         pot = pot - charge(3)*log(d1)
      endif
C
      if (ifdouble.eq.1) then
         d2 = d1/rr
         u11 = d1*(-1.0d0+d2*r1**2)
         u12 = d1*d2*r1*r2
         u22 = d1*(-1.0d0+d2*r2**2)
         u33 = -1.0d0/rr
         pot2 = dipstr(1)*dipvec(1)*u11 +
     1          dipstr(2)*dipvec(2)*u22 +
     1          dipstr(3)*dipvec(3)*u33 +
     1          (dipstr(1)*dipvec(2)+dipstr(2)*dipvec(1))*u12
         rnudotf = dipstr(1)*dipvec(1)+
     1             dipstr(2)*dipvec(2)+
     1             dipstr(3)*dipvec(3)
         pot = pot + 2*(rlame(2)*pot2 + rnudotf*u33*rlame(1))
      endif
C
      pot = pot*scfac
      return
      end
C
C
C
C***********************************************************************
      subroutine l3dformtaintkercf90(ier,scfac,source,ifcharge,charge,
     1           ifdouble,dipstr,dipvec,
     1           ns,rlame,center,radius,nterms,local,local2)
C***********************************************************************
C
C     This subrotine forms the local expansion due to a collection 
C     of Mindlin C sources.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     ier    = error return code for memory allocation inside routine
C     scfac     = scale parameter for local expansion.
C     source    = source locations
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of first C   s.h. expansion
C     local2 = coefficients of second C  s.h. expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ier,lw,lused,ns,nquad,nquadm,nterms
      integer ifcharge,ifdouble
      integer iphival,ifinit,iprojloc,iused,iwts,ixnodes
      integer lynm,iynm,lphival
      real *8 source(3,ns),center(3)
      real *8 scfac,radius,rlame(2)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 local2(0:nterms,-nterms:nterms)
      real *8, allocatable :: w(:)
c
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
      call l3devalsphereintkerc1(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,rlame,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call projloc3d(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local,w(iprojloc),w(iynm))
c
      call l3drescale2(nterms,local,radius,scfac)
c
      call l3devalsphereintkerc2(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,rlame,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call projloc3d(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local2,w(iprojloc),w(iynm))
c
      call l3drescale2(nterms,local2,radius,scfac)
      return
      end
c
c
c
C***********************************************************************
      subroutine l3devalsphereintkerc1(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,rlame,
     1     phival,center,radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine tabulates the Mindlin C1 harmonic function 
C     due to a collection of sources/charges
C     at quadrature nodes on a sphere of radius RADIUS,
C     centered at CENTER.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source    = source locations
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C     nquad     = number of quadrature nodes in theta direction
C     nquadm    = number of quadrature nodes in phi direction
C        [total number of nodes on target sphere is nquad*nquadm]
C     xnodes    = Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C---------------------------------------------------------------------
      implicit none
      integer i,jj,kk,ns,nterms,nquad,nquadm
      integer ifcharge,ifdouble
      real *8 source(3,ns),targ(3), center(3)
      real *8 xnodes(nquad)
      real *8 rlame(2)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      real *8 radius,pi,ctheta,stheta,phi,cosphi,sinphi
      complex *16 phival(nquad,nquadm)
      complex *16 pot,fld(3)
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
            call intkerc1(rlame,source(1,i),ifcharge,charge(1,i),
     1           ifdouble,dipstr(1,i),dipvec(1,i),targ,pot)
            phival(jj,kk) = phival(jj,kk) + pot
         enddo
      enddo
      enddo
      return
      end
C
C***********************************************************************
      subroutine intkerc1(rlame,source,ifcharge,charge,
     1           ifdouble,dipstr,dipvec,targ,pot)
C***********************************************************************
C
C     This subroutine computes the harmonic Mindlin C1 kernel
C     at a target point.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     rlame     = Lame coefficients
C     source    = source location
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     targ      = target location
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     pot       =  value of Mindlin B pre-harmonic potential 
C---------------------------------------------------------------------

      implicit none
      integer ifcharge,ifdouble
      integer nquad,nquadm,jj,kk
      real *8 source(3),targ(3)
      real *8 rlame(2)
      real *8 charge(3)
      real *8 dipstr(3)
      real *8 dipvec(3)
      real *8 alpha,r1,r2,r3,rr,rr3,rr5
      real *8 v,v1,v2,v3,scfac,d3,d5,v11,v12,v13,v22,v23,v33
      real *8 pot2,rfac1,rnudotf
      complex *16 pot
C
      pot = 0.0d0
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      r1 = targ(1)-source(1)
      r2 = targ(2)-source(2)
      r3 = targ(3)-source(3)
      rr = sqrt(r1*r1+r2*r2+r3*r3)
      rr3 = rr*rr*rr
      v = 1.0d0/rr
      v1 = -r1/rr3
      v2 = -r2/rr3
      v3 = -r3/rr3
      if (ifcharge.eq.1) then
         scfac  = alpha*source(3)
         pot = pot + scfac*(charge(1)*v1+charge(2)*v2-charge(3)*v3)
         pot = pot - (2.0d0-alpha)*charge(3)*v
      endif
c
      if (ifdouble.eq.1) then
         rr5 = rr3*rr*rr
         d3 = 1.0d0/rr3
         d5 = 1.0d0/rr5
         v11 = -d3 + 3*d5*r1**2
         v12 = 3*d5*r1*r2
         v13 = 3*d5*r1*r3
         v22 = -d3 + 3*d5*r2**2
         v23 = 3*d5*r2*r3
         v33 = -d3 + 3*d5*r3**2
         rfac1 = 2.0d0*alpha*source(3)*rlame(2)
         pot2 = -rfac1*(dipstr(1)*dipvec(1)*v11 +
     1                 dipstr(2)*dipvec(2)*v22 +
     1                 dipstr(3)*dipvec(3)*v33 +
     1          (dipstr(1)*dipvec(2)+dipstr(2)*dipvec(1))*v12 -
     1          (dipstr(1)*dipvec(3)+dipstr(3)*dipvec(1))*v13 -
     1          (dipstr(2)*dipvec(3)+dipstr(3)*dipvec(2))*v23)
         rfac1 = (2.0d0-2.0d0*alpha)*rlame(2)
         rnudotf = dipstr(1)*dipvec(1)+
     1             dipstr(2)*dipvec(2)+
     1             dipstr(3)*dipvec(3)
         pot2 = pot2 +rfac1*(
     1          (dipstr(1)*dipvec(3)+dipstr(3)*dipvec(1))*v1 +
     1          (dipstr(2)*dipvec(3)+dipstr(3)*dipvec(2))*v2 -
     1          (dipstr(3)*dipvec(3)+dipstr(3)*dipvec(3))*v3)
         pot2 = pot2 +rlame(1)*2.0d0*(alpha-1)*rnudotf*v3
      endif
      pot = pot + pot2
c
      return
      end
c
c
C***********************************************************************
      subroutine l3devalsphereintkerc2(source,ifcharge,charge,ifdouble,
     1     dipstr,dipvec,ns,rlame,
     1     phival,center,radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine tabulates the Mindlin C2 harmonic function 
C     due to a collection of sources/charges
C     at quadrature nodes on a sphere of radius RADIUS,
C     centered at CENTER.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source    = source locations
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C     nquad     = number of quadrature nodes in theta direction
C     nquadm    = number of quadrature nodes in phi direction
C        [total number of nodes on target sphere is nquad*nquadm]
C     xnodes    = Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C---------------------------------------------------------------------
      implicit none
      integer i,jj,kk,ns,nterms,nquad,nquadm
      integer ifcharge,ifdouble
      real *8 source(3,ns),targ(3), center(3)
      real *8 xnodes(nquad)
      real *8 rlame(2)
      real *8 charge(3,ns)
      real *8 dipstr(3,ns)
      real *8 dipvec(3,ns)
      real *8 radius,pi,ctheta,stheta,phi,cosphi,sinphi
      complex *16 phival(nquad,nquadm)
      complex *16 pot,fld(3)
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
            call intkerc2(rlame,source(1,i),ifcharge,charge(1,i),
     1           ifdouble,dipstr(1,i),dipvec(1,i),targ,pot)
            phival(jj,kk) = phival(jj,kk) + pot
         enddo
      enddo
      enddo
      return
      end
C
C***********************************************************************
      subroutine intkerc2(rlame,source,ifcharge,charge,
     1           ifdouble,dipstr,dipvec,targ,pot)
C***********************************************************************
C
C     This subroutine computes the harmonic Mindlin C2 kernel
C     at a target point.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     rlame     = Lame coefficients
C     source    = source location
C     ifcharge  = flag for single layer source
C     charge    = single layer force vector
C     ifdouble  = flag for double layer source
C     dipstr    = double layer force vector
C     dipvec    = double layer normal vector
C     targ      = target location
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     pot       =  value of Mindlin B pre-harmonic potential 
C---------------------------------------------------------------------

      implicit none
      integer nquad,nquadm,jj,kk
      integer ifcharge,ifdouble
      real *8 source(3),targ(3)
      real *8 rlame(2)
      real *8 charge(3)
      real *8 dipstr(3)
      real *8 dipvec(3)
      real *8 alpha,r1,r2,r3,rr,rr3,rr5
      real *8 d3,d5,v,v1,v2,v3,v11,v12,v13,v22,v23,v33,rfac1,pot2
      complex *16 pot
C
      pot = 0.0d0
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      r1 = targ(1)-source(1)
      r2 = targ(2)-source(2)
      r3 = targ(3)-source(3)
      rr = sqrt(r1*r1+r2*r2+r3*r3)
      rr3 = rr*rr*rr
      v = 1.0d0/rr
      v1 = -r1/rr3
      v2 = -r2/rr3
      v3 = -r3/rr3
      if (ifcharge.eq.1) then
         pot = pot + (2.0d0-alpha)*
     1         (charge(1)*v1+charge(2)*v2-charge(3)*v3)
      endif
c
      if (ifdouble.eq.1) then
         rr5 = rr3*rr*rr
         d3 = 1.0d0/rr3
         d5 = 1.0d0/rr5
         v11 = -d3 + 3*d5*r1**2
         v12 = 3*d5*r1*r2
         v13 = 3*d5*r1*r3
         v22 = -d3 + 3*d5*r2**2
         v23 = 3*d5*r2*r3
         v33 = -d3 + 3*d5*r3**2
         rfac1 = 2.0d0*(2.0d0-alpha)*rlame(2)
         pot2 = rfac1*(dipstr(1)*dipvec(1)*v11 +
     1          dipstr(2)*dipvec(2)*v22 +
     1          dipstr(3)*dipvec(3)*v33 +
     1          (dipstr(1)*dipvec(2)+dipstr(2)*dipvec(1))*v12 -
     1          (dipstr(1)*dipvec(3)+dipstr(3)*dipvec(1))*v13 -
     1          (dipstr(2)*dipvec(3)+dipstr(3)*dipvec(2))*v23)
      endif
      pot = pot - pot2
      return
      end
c
c
c**********************************************************************
      subroutine mindlintaevalb(rscale,center,local,nterms,
     1		ztarg,fld,ifhess,gradfld)
c**********************************************************************
c     This subroutine evaluates the displacement (fld) and gradient
c     of displacement (gradfld) of the 
c     local expansion due to Mindlin B sources.
c
c     fld = -gradient
c     gradfld = Hessian with gradfld(i,j) = dx_i( fld(j)).
c
c     where rscale defines scaling parameter.     
c
c     Subroutine for computing the displacement and gradient of 
c     displacement for a Mindlin B local expansion.
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale =    scaling parameter 
c     center =    expansion center
c     local  =    local expansion in 2d matrix format
c     nterms =    order of the multipole expansion
c     ztarg  =    target location
C     ifhess    =  flag to compute gradient of displacement
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     fld       =    -gradient at ztarg
c     gradfld   =    gradient of fld (if requested)
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ifhess,ier,iffld
      real *8 center(3),ztarg(3)
      real *8 rscale
      complex *16 pot,fld(3)
      complex *16 hess(6)
      complex *16 gradfld(3,3)
      complex *16 local(0:nterms,-nterms:nterms)
c
      iffld = 1
      call l3dtaevalhess(rscale,center,local,nterms,ztarg,pot,
     1     iffld,fld,ifhess,hess,ier)
      fld(1) = -fld(1)
      fld(2) = -fld(2)
c
      if (ifhess.eq.1) then
         gradfld(1,1) = -hess(1)
         gradfld(2,2) = -hess(2)
         gradfld(3,3) = hess(3)
         gradfld(1,2) = -hess(4)
         gradfld(2,1) = -hess(4)
         gradfld(1,3) = hess(5)
         gradfld(3,1) = -hess(5)
         gradfld(2,3) = hess(6)
         gradfld(3,2) = -hess(6)
      endif
      return
      end
c
c**********************************************************************
      subroutine mindlintaevalc(rscale,center,local,local2,nterms,
     1		ztarg,fld,ifhess,gradfld)
c**********************************************************************
c     This subroutine evaluates the displacement (fld) and gradient
c     of displacement (gradfld) of the 
c     local expansion due to Mindlin C sources.
c
c     gradfld = Hessian with gradfld(i,j) = dx_i( fld(j)).
c
c     Subroutine for computing the displacement and gradient of 
c     displacement for a Mindlin B local expansion.
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale =    scaling parameter for local expansion
c     center =    expansion center
c     local  =    first local expansion in 2d matrix format
c     local2  =    second local expansion in 2d matrix format
c     nterms =    order of the multipole expansion
c     ztarg  =    target location
C     ifhess    =  flag to compute gradient of displacement
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     fld       =    -gradient at ztarg
c     gradfld   =    gradient of fld (if requested)
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ifhess,ier,iffld
      real *8 center(3),ztarg(3)
      real *8 rscale
      complex *16 pot,fld(3),fld2(3)
      complex *16 hess(6)
      complex *16 gradfld(3,3)
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 local2(0:nterms,-nterms:nterms)
c
c     
      iffld = 1
      call l3dtaevalhess(rscale,center,local,nterms,ztarg,pot,
     1     iffld,fld,ifhess,hess,ier)
      fld(1) = -fld(1)
      fld(2) = -fld(2)
      fld(3) = -fld(3)
      call l3dtaeval(rscale,center,local2,nterms,ztarg,pot,
     1     ifhess,fld2,ier)
      fld(3) = fld(3) - pot
c
      if (ifhess.eq.1) then
         gradfld(1,1) = -hess(1)
         gradfld(2,2) = -hess(2)
         gradfld(3,3) = -hess(3) - fld2(3)
         gradfld(1,2) = -hess(4)
         gradfld(2,1) = -hess(4)
         gradfld(1,3) = -hess(5) - fld2(1)
         gradfld(3,1) = -hess(5)
         gradfld(2,3) = -hess(6) - fld2(2)
         gradfld(3,2) = -hess(6)
      endif
      return
      end
C
C
C
c**********************************************************************
      subroutine mindlintaevalb_trunc(rscale,center,local,nterms,
     1		ztarg,fld,ifhess,gradfld,
     $          scarray_local,wlege,nlege)
c**********************************************************************
c     This subroutine evaluates the displacement (fld) and gradient
c     of displacement (gradfld) of the 
c     local expansion due to Mindlin B sources.
c
c     fld = -gradient
c     gradfld = Hessian with gradfld(i,j) = dx_i( fld(j)).
c
c     where rscale defines scaling parameter.     
c
c     Subroutine for computing the displacement and gradient of 
c     displacement for a Mindlin B local expansion.
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale =    scaling parameter 
c     center =    expansion center
c     local  =    local expansion in 2d matrix format
c     nterms =    order of the multipole expansion
c     ztarg  =    target location
C     ifhess    =  flag to compute gradient of displacement
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     fld       =    -gradient at ztarg
c     gradfld   =    gradient of fld (if requested)
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ifhess,ier,iffld
      real *8 center(3),ztarg(3)
      real *8 rscale
      complex *16 pot,fld(3)
      complex *16 hess(6)
      complex *16 gradfld(3,3)
      complex *16 local(0:nterms,-nterms:nterms)
      integer nlege
      real *8 scarray_local(*),wlege(*)
c
      iffld = 1
c      call l3dtaevalhess(rscale,center,local,nterms,ztarg,pot,
c     1     iffld,fld,ifhess,hess,ier)
      if( ifhess .eq. 1 ) then
      call l3dtaevalhessd_trunc(rscale,center,local,nterms,ztarg,pot,
     1     iffld,fld,ifhess,hess,scarray_local,wlege,nlege)
      else
      call l3dtaeval_trunc(rscale,center,local,nterms,nterms,
     $     ztarg,pot,iffld,fld,wlege,nlege,ier)
      endif
      fld(1) = -fld(1)
      fld(2) = -fld(2)
c
      if (ifhess.eq.1) then
         gradfld(1,1) = -hess(1)
         gradfld(2,2) = -hess(2)
         gradfld(3,3) = hess(3)
         gradfld(1,2) = -hess(4)
         gradfld(2,1) = -hess(4)
         gradfld(1,3) = hess(5)
         gradfld(3,1) = -hess(5)
         gradfld(2,3) = hess(6)
         gradfld(3,2) = -hess(6)
      endif
      return
      end
c
c**********************************************************************
      subroutine mindlintaevalc_trunc
     $     (rscale,center,local,local2,nterms,
     1		ztarg,fld,ifhess,gradfld,
     $          scarray_local,wlege,nlege)
c**********************************************************************
c     This subroutine evaluates the displacement (fld) and gradient
c     of displacement (gradfld) of the 
c     local expansion due to Mindlin C sources.
c
c     gradfld = Hessian with gradfld(i,j) = dx_i( fld(j)).
c
c     Subroutine for computing the displacement and gradient of 
c     displacement for a Mindlin B local expansion.
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale =    scaling parameter for local expansion
c     center =    expansion center
c     local  =    first local expansion in 2d matrix format
c     local2  =    second local expansion in 2d matrix format
c     nterms =    order of the multipole expansion
c     ztarg  =    target location
C     ifhess    =  flag to compute gradient of displacement
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     fld       =    -gradient at ztarg
c     gradfld   =    gradient of fld (if requested)
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ifhess,ier,iffld
      real *8 center(3),ztarg(3)
      real *8 rscale
      complex *16 pot,fld(3),fld2(3)
      complex *16 hess(6)
      complex *16 gradfld(3,3)
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 local2(0:nterms,-nterms:nterms)
c
      integer nlege
      real *8 scarray_local(*),wlege(*)
c
c     
      iffld = 1
c      call l3dtaevalhess(rscale,center,local,nterms,ztarg,pot,
c     1     iffld,fld,ifhess,hess,ier)
      if(ifhess .eq. 1 ) then
      call l3dtaevalhessd_trunc(rscale,center,local,nterms,ztarg,pot,
     1     iffld,fld,ifhess,hess,scarray_local,wlege,nlege)
      else
      call l3dtaeval_trunc(rscale,center,local,nterms,nterms,
     $     ztarg,pot,iffld,fld,wlege,nlege,ier)
      endif
      fld(1) = -fld(1)
      fld(2) = -fld(2)
      fld(3) = -fld(3)
c      call l3dtaeval
c     $   (rscale,center,local2,nterms,ztarg,pot,
c     1     ifhess,fld2,ier)
      call l3dtaeval_trunc
     $   (rscale,center,local2,nterms,nterms,ztarg,pot,
     1     ifhess,fld2,wlege,nlege,ier)
      fld(3) = fld(3) - pot
c
      if (ifhess.eq.1) then
         gradfld(1,1) = -hess(1)
         gradfld(2,2) = -hess(2)
         gradfld(3,3) = -hess(3) - fld2(3)
         gradfld(1,2) = -hess(4)
         gradfld(2,1) = -hess(4)
         gradfld(1,3) = -hess(5) - fld2(1)
         gradfld(3,1) = -hess(5)
         gradfld(2,3) = -hess(6) - fld2(2)
         gradfld(3,2) = -hess(6)
      endif
      return
      end
C
C
C
C***********************************************************************
      subroutine l3dformtaintkerbf90sc(ier,scfac,rscale,source,charge1,
     1           ns,rlame,center,radius,nterms,local)
C***********************************************************************
C
C     This subroutine forms the local expansion due to a collection 
C     of Mindlin B sources.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     ier    = error return code for memory allocation inside routine
C     scfac     = scale parameter for local expansion.
C     rscale    = scaling parameter for original multipole expansion
C     source    = source locations
C     charge1   = complex pt source strengths
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of s.h. expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ier,lw,lused,ns,nquad,nquadm,nterms
      integer ifcharge,ifdouble
      integer ifinit
      integer iphival,lphival,ixnodes,iwts,iynm,lynm,iprojloc,iused
      real *8 source(3,ns),center(3)
      real *8 rscale,scfac,radius,rlame(2)
      complex *16 charge1(ns)
      complex *16 local(0:nterms,-nterms:nterms)
      real *8, allocatable :: w(:)
c
c     carve up workspace and compute allocation needed
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
      call l3devalsphereintkerbsc(source,charge1,ns,rlame,rscale,
     1     w(iphival),center,radius,nterms,nquad,nquadm,w(ixnodes))
c
      call projloc3d(nterms,nterms,nquad,w(ixnodes),w(iwts),
     1     w(iphival),local,w(iprojloc),w(iynm))
c
      call l3drescale2(nterms,local,radius,scfac)
      return
      end
c
c
c
C***********************************************************************
      subroutine l3devalsphereintkerbsc(source,charge1,ns,rlame,
     1     rscale,phival,center,radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine tabulates a harmonic function defined by 
C     a user-provided subroutine and a collection of sources/charges
C     at quadrature nodes on a sphere of radius RADIUS,
C     centered at CENTER.
C
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source    = source locations
C     charge1    = single layer force vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     center    = expansion center
C     radius    = sphere radius
C     nterms    = order of spherical harmonic expansion
C     nquad     = number of quadrature nodes in theta direction
C     nquadm    = number of quadrature nodes in phi direction
C        [total number of nodes on target sphere is nquad*nquadm]
C     xnodes    = Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit none
      integer nterms,ns,ifcharge,ifdouble
      integer nquad,nquadm,jj,kk,i
      real *8 source(3,ns),targ(3), center(3)
      real *8 xnodes(nquad)
      real *8 rlame(2),radius
      real *8 rscale,r1,r2,r3,rr,rkerb,scfac,alpha
      real *8 cosphi,ctheta,phi,pi,sinphi,stheta
      complex *16 charge1(ns)
      complex *16 phival(nquad,nquadm)
      complex *16 pot
C
      pi = 4.0d0*datan(1.0d0)
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      scfac  = (1.0d0-alpha)/alpha
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
         do i = 1,ns
            r1 = targ(1)-source(1,i)
            r2 = targ(2)-source(2,i)
            r3 = targ(3)-source(3,i)
            rr = sqrt(r1*r1+r2*r2+r3*r3)
            rkerb = -r3*dlog(rr-r3) - rr
            pot = charge1(i)*rkerb
            phival(jj,kk) = phival(jj,kk) + pot
         enddo
         phival(jj,kk) = phival(jj,kk)*scfac*rscale
      enddo
      enddo
      return
      end
C
