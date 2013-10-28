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
c    $Date: 2010-05-10 12:12:17 -0400 (Mon, 10 May 2010) $
c    $Revision: 958 $
c
c
c    Subroutine for converting multipole expansion to corresponding
c    set of point charges on sphere of given radius.
c
C***********************************************************************
      subroutine mptoslpborder(mpole,nterms,center,rscale,radius,
     1           nquadm,source,charge)
C***********************************************************************
C     This subroutine computes a set of charges on the sphere
C     that match the <<border>> of the multipole expansion:
C     that is the terms of the form 
C     M(n,n), M(n+1,n) for n = 0,nterms.
C     There are only 2*nquadm = 0(nterms) such discrete charges required.
C     In the Mindlin B multipole expansion, the difference can then
C     be handled analytically, by using the dz^{-2}
C     antiderivative.
C
C     Used in debugging faster mpeval, where we convert mpole
C     to dz^{-2} mpole for bulk of expansion but border needs to 
C     be handled separately.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     mpole  = multipole expansion
C     nterms = order of spherical harmonic expansion
C     center = center of expansion
C     rscale = multipole scaling parameter
C     radius = radius of sphere for charge discretization
C     nquadm =  number of quadrature nodes 
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           source(3,2*nquadm) = charge locations 
C           charge(2*nquadm) = complex charge values
C
C***********************************************************************
      implicit none
      integer nterms,nquadm
      integer l,m,j,k
      real *8 center(3),pi,ctheta,stheta,rfac,rmul,rscale,radius
      real *8 cosphi,sinphi,phi,h
      real *8 source(3,2*nquadm)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 chat2(-100:100),chat1(-100:100),zf1,zf2
      complex *16 ephi,imag,emul,sum,zmul,zfac
      complex *16 charge(2*nquadm)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
      ctheta = 1.0d0/sqrt(3.0d0)
      stheta = dsqrt(2.0d0/3.0d0)
      call ylgndr(nterms,ctheta,ynm)
      rfac = 1.0d0/(radius*rscale)
      rmul = 1.0d0
      do j = 0,nterms
ccc         zf1 = mpole(j,j)/(ynm(j,j)*stheta*(radius**(j)))
ccc         zf1 = zf1*sqrt(2*j+1.0d0)/(rscale**j)
         zf1 = mpole(j,j)/(ynm(j,j)*stheta)*rmul
         rmul = rmul*rfac
         zf1 = zf1*sqrt(2*j+1.0d0)
         if (j+1.le. nterms) then
ccc            zf2 = mpole(j+1,j)/(ynm(j+1,j)*stheta*(radius**(j+1)))
ccc            zf2 = zf2*sqrt(2*(j+1)+1.0d0)/(rscale**(j+1))
            zf2 = mpole(j+1,j)/(ynm(j+1,j)*stheta)*rmul
            zf2 = zf2*sqrt(2*(j+1)+1.0d0)
         else
            zf2 = 0.0d0
         endif
         chat1(j) = (zf1+zf2)/2.0d0
         chat2(j) = (zf1-zf2)/2.0d0
      enddo
      do j = 1,nterms
         chat1(-j) = dconjg(chat1(j))
         chat2(-j) = dconjg(chat2(j))
      enddo
c
      h = 2*pi/nquadm
      do j = 1,nquadm
	 phi = (j-1)*h
         charge(j) = 0.0d0
         charge(j+nquadm) = 0.0d0
         zmul = cdexp(imag*(-nterms-1)*phi)
         zfac = cdexp(imag*phi)
         do k = -nterms,nterms
            zmul = zmul*zfac
ccc            charge(j) = charge(j) + chat1(k)*cdexp(imag*k*phi)
            charge(j) = charge(j) + chat1(k)*zmul
            charge(j+nquadm) = charge(j+nquadm) + 
ccc     $                  chat2(k)*cdexp(imag*k*phi)
     $                  chat2(k)*zmul
         enddo
         charge(j) = charge(j)*stheta/nquadm
         charge(j+nquadm) = charge(j+nquadm)*stheta/nquadm
      enddo
      do j = 1,nquadm
	 phi = (j-1)*h
	 cosphi = dcos(phi)
	 sinphi = dsin(phi)
	 source(1,j) = center(1)+radius*stheta*cosphi
	 source(2,j) = center(2)+radius*stheta*sinphi
	 source(3,j) = center(3)+radius*ctheta
         source(1,j+nquadm) = source(1,j)
         source(2,j+nquadm) = source(2,j)
         source(3,j+nquadm) = center(3)-radius*ctheta
      enddo
c
      return
      end
c
c
