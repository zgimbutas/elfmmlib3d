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
C***********************************************************************
      subroutine projloc3d(nterms,ldl,nquad,xnodes,wts,phival,local,
     1           marray,ynm)
C***********************************************************************
C
C     Usage:
C
C           compute spherical harmonic expansion on unit sphere
C           of function tabulated at nquad*nquad grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl    = dimension parameter for local expansion
C           nquad  = number of quadrature nodes in each direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           marray  = workspace of dimension (nquad,-nterms:nterms)
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C
C    NOTE:
C
C    yrecursion.f produces Ynm with a nonstandard scaling:
C    (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C    is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C   In the first loop below, you see
C
Cccc	    marray(jj,m) = sum*2*pi/nquad
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C***********************************************************************
      implicit none
      integer nterms,ldl,nquad,ier
      integer l,m,jj,kk
      real *8 wts(nquad),xnodes(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 cthetaj,pi
      complex *16 zk,phival(nquad,nquad)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 marray(nquad,-nterms:nterms)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
         enddo
      enddo
c
c     create marray (intermediate array)
c
      emul1 = cdexp(imag*2*pi/nquad)
      emul = cdexp(-imag*2*nterms*pi/nquad)
      do m=-nterms,nterms
cc	    emul = cdexp(imag*m*2*pi/nquad)
	 do jj=1,nquad
	    sum = 0
	    ephi = emul
	    do kk = 1,nquad
               sum = sum + phival(jj,kk)*dconjg(ephi)
	       ephi = ephi*emul
            enddo
ccc	    marray(jj,m) = sum*2*pi/nquad
	    marray(jj,m) = sum/(2*nquad)
         enddo
	 emul = emul*emul1
      enddo
c
c     get local exp
c
      do jj=1,nquad
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-nterms,nterms
	    zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
