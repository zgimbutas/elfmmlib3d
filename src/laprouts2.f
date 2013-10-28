cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
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
c      This file contains the basic subroutines for 
c      forming and evaluating multipole expansions.
c
c      Remarks on scaling conventions.
c
c      1)  Far field and local expansions are consistently rscaled as
c              
c
c          M_n^m (scaled) = M_n^m / rscale^(n)  so that upon evaluation
c
c          the field is  sum   M_n^m (scaled) * rscale^(n) / r^{n+1}.
c
c          L_n^m (scaled) = L_n^m * rscale^(n)  so that upon evaluation
c
c          the field is  sum   L_n^m (scaled) / rscale^(n) * r^{n}.
c
c
c      2) There are many definitions of the spherical harmonics,
c         which differ in terms of normalization constants. We
c         adopt the following convention:
c
c         For m>0, we define Y_n^m according to 
c
c         Y_n^m = \sqrt{2n+1} \sqrt{\frac{ (n-m)!}{(n+m)!}} \cdot
c                 P_n^m(\cos \theta)  e^{i m phi} 
c         and
c 
c         Y_n^-m = dconjg( Y_n^m )
c    
c         We omit the Condon-Shortley phase factor (-1)^m in the 
c         definition of Y_n^m for m<0. (This is standard in several
c         communities.)
c
c         We also omit the factor \sqrt{\frac{1}{4 \pi}}, so that
c         the Y_n^m are orthogonal on the unit sphere but not 
c         orthonormal.  (This is also standard in several communities.)
c         More precisely, 
c
c                 \int_S Y_n^m Y_n^m d\Omega = 4 \pi. 
c
c         Using our standard definition, the addition theorem takes 
c         the simple form 
c
c         1/r = 
c         \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)/ (|T|^(n+1)) 
c
c         1/r = 
c         \sum_n \sum_m  |S|^n  Ylm*(S)    Ylm(T)     / (|T|^(n+1)) 
c                               -------    ------
c                               sqrt(2n+1) sqrt(2n+1)
c
c        In the Laplace library (this library), we incorporate the
c        sqrt(2n+1) factor in both forming and evaluating multipole
c        expansions.
c
c-----------------------------------------------------------------------
c
c      f90 version, using allocate
c
c
c      L3DMPEVAL: computes potential and -grad(potential)
c                 due to a multipole expansion.
c                 (calls L3DMPEVAL0)
c
c      L3DFORMMP: creates multipole expansion (outgoing) due to 
c                 a collection of sources.
c                 (calls FORMMPONEL3D/FORMMPL3D0 for each source)
c
c      CART2POLARL: utility function.
c                  converts Cartesian coordinates into polar
c                  representation needed by other routines.
c
c      L3DTAEVAL: computes potential and -grad(potential) 
c                  due to local expansion.
c                 (calls TAEVALL3D0)
c
c      L3DFORMTA: creates local expansion due to 
c                 a collection of sources.
c                 (calls FORMTAL3DONE/FORMTAL3D0 for each source)
c
c      LPOTFLD3DALL:  direct calculation for a set of charge sources
c      LPOTFLD3D : direct calculation for a single charge source
c
c      L3DADD   : adds one expansion to another
c
c      L3DFORMMP_DP: creates multipole expansion (outgoing) due to 
c                 a collection of dipoles.
c                 (calls FORMMPONEL3D_DP/FORMMPL3D0_DP for each source)
c
c      LPOTFLD3DALL_DP:  direct calculation for a set of dipole sources
c      LPOTFLD3D_DP : direct calculation for a single dipole source
c
c      L3DFORMTA_DP: creates local expansion (incoming) due to 
c                 a collection of dipoles.
c                 (calls L3DFORMTA1_DP and L3DFORMTA0_DP)
c
c
c
c
c
c**********************************************************************
      subroutine l3dmpeval(rscale,center,mpole,nterms,ztarg,
     1		pot,iffld,fld,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi) / r^{n+1} / sqrt(2n+1)
c             n   m
c
c     fld = -gradient(pot) if iffld = 0.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    -gradient at ztarg (if requested)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c-----------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ier
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lused
      real *8 rscale
      real *8 center(3),ztarg(3)
      real *8, allocatable :: w(:)
      complex *16 pot,fld(3)
      complex *16 mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
c     for Ynm and Ynm'
c
      lpp=(nterms+1)**2+5
      ipp=1
      ippd = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+3)+5 
c
      ifr=iephi+lephi
      ifrder=ifr+(nterms+3)
      lused=ifrder+(nterms+3)
      allocate(w(lused))
c
      call l3dmpeval0(rscale,center,mpole,nterms,ztarg,
     1	   pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifr),w(ifrder))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dmpeval0(rscale,center,mpole,nterms,
     1		ztarg,pot,iffld,fld,ynm,ynmd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dmpeval for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,iffld
      integer i,l,n,m
      real *8 rscale
      real *8 center(3),ztarg(3),zdiff(3)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 fr(0:nterms+1)
      real *8 frder(0:nterms+1)
      complex *16 pot,fld(3),ephi1,ur,utheta,uphi,ux,uy,uz
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
c
      real *8 phi,cphi,sphi,theta,ctheta,stheta,thetax,thetay,thetaz
      real *8 rs,rx,ry,rz,d,r,phix,phiy,phiz,done
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = 1.0d0/r
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=done
      ephi(1)=ephi1
      cphi = dreal(ephi1)
      sphi = dimag(ephi1)
      ephi(-1)=dconjg(ephi1)
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      do i=0,nterms+1
         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
      if (iffld.eq.1) then
         rx = stheta*cphi
         thetax = ctheta*cphi/r
         phix = -sphi/r
         ry = stheta*sphi
         thetay = ctheta*sphi/r
         phiy = cphi/r
         rz = ctheta
         thetaz = -stheta/r
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions
c     and scale by 1/sqrt(2l+1)
c
      if (iffld.eq.1) then
         call ylgndr2s(nterms,ctheta,ynm,ynmd)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               ynm(l,m) = ynm(l,m)*rs
               ynmd(l,m) = ynmd(l,m)*rs
            enddo
         enddo
      else        
         call ylgndr(nterms,ctheta,ynm)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               ynm(l,m) = ynm(l,m)*rs
            enddo
         enddo
      endif
c
c     initialize computed values.
c
      if (iffld.eq.1) then
         ur = mpole(0,0)*frder(0)
         utheta = 0.0d0
         uphi = 0.0d0
      endif
      pot=mpole(0,0)*fr(0)
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms
	    pot=pot+mpole(n,0)*fr(n)*ynm(n,0)
	    ur = ur + frder(n)*ynm(n,0)*mpole(n,0)
	    utheta = utheta -mpole(n,0)*fr(n)*ynmd(n,0)*stheta
	    do m=1,n
	       ztmp1=fr(n)*ynm(n,m)*stheta
	       ztmp2 = mpole(n,m)*ephi(m) 
	       ztmp3 = mpole(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + frder(n)*ynm(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fr(n)*ynmd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fr(n)*ynm(n,m)*ztmpsum
            enddo
         enddo
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms
	    pot=pot+mpole(n,0)*fr(n)*ynm(n,0)
	    do m=1,n
	       ztmp1=fr(n)*ynm(n,m)
	       ztmp2 = mpole(n,m)*ephi(m) + mpole(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine l3dformmp(ier,rscale,sources,charge,ns,center,
     1                  nterms,mpole)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		                   ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c    
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, ier1, lused
      real *8 center(3),sources(3,ns)
      real *8 rscale,rs
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,charge(ns)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      ier = 0
      do i = 1, ns
         call l3dformmp1(ier1,rscale,sources(1,i),charge(i),center,
     1        nterms,mpole)
      enddo
      if (ier1.ne.0) ier = ier1
c
c     scale by 1/sqrt(2l+1)
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            mpole(l,m) = mpole(l,m)*rs
         enddo
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1(ier,rscale,source,charge,center,
     1		nterms,mpole)
c**********************************************************************
c
c     This subroutine creates the h-expansion about CENTER
c     due to a charge located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     charge  : complex charge strength
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c                            
c     mpole   : coeffs of the h-expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,ippd,iephi,lephi,ifrder,lfrder,ifr,lused
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 charge
c
c     compute work space components:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifrder=iephi+lephi
      lfrder=2*(nterms+3)
c
      ifr=ifrder+lfrder
      lused=ifr+lfrder

      allocate(w(lused))
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp0(rscale,source,charge,center,nterms,
     1		mpole,w(ipp),w(ippd),w(iephi),w(ifr),
     2          w(ifrder))
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0(rscale,source,charge,center,
     1		nterms,mpole,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dformmp1 for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,n,m
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      real *8 r,theta,phi,d,ctheta,stheta,cphi,sphi,dtmp
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 charge
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16  fr(0:nterms+1),frder(0:nterms+1)
      complex *16  ztmp,z
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = r
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      d = d*rscale
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i)=i*fr(i-1)
      enddo
c
c     get the associated Legendre functions:
c
      call ylgndr(nterms,ctheta,pp)
ccc      call ylgndr2s(nterms,ctheta,pp,ppd)
ccc      call prinf(' after ylgndr with nterms = *',nterms,1)
ccc      call prinm2(pp,nterms)
c
c     multiply all fr's by charge strength.
c
      do n = 0,nterms
         fr(n) = fr(n)*charge
      enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r =  
c          \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)  / (|T|)^{n+1}
c
c     so contribution is |S|^n times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      mpole(0,0)= mpole(0,0) + fr(0)
      do n=1,nterms
         dtmp=pp(n,0)
         mpole(n,0)= mpole(n,0) + dtmp*fr(n)
         do m=1,n
            ztmp=pp(n,m)*fr(n)
            mpole(n, m)= mpole(n, m) + ztmp*dconjg(ephi(m))
            mpole(n,-m)= mpole(n,-m) + ztmp*dconjg(ephi(-m))
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dtaeval(rscale,center,locexp,nterms,
     1		ztarg,pot,iffld,fld,ier)
c**********************************************************************
c
c     This subroutine evaluates a local expansion centered at CENTER
c     at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) r^n Y_nm(theta,phi) / sqrt(2n+1)
c             n   m
c
c     The reason for including the sqrt(2n+1) scaling has to do with
c     the addition theorem for 1/r. The term 1/(2n+1) is needed
c     and we put half the weight on the local evaluation and half the 
c     weight on the expansion formation...
c
c     The addition theorem for exp(ikr)/r does not require the 
c     1/(2n+1) scaling - it appears in the definitions of the 
c     Bessel and Hankel functions.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale     : scaling parameter used in forming expansion
c                                   (see l3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the j-expansion
c     nterms     : order of the h-expansion
c     ztarg      : target vector
c     iffld      : flag for gradient computation
c		        iffld=0  - gradient is not computed
c		        iffld=1  - gradient is computed
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot        : potential at ztarg(3)
c     fld        : -gradient at ztarg (if requested)
c     lused      : amount of work space "w" used
c     ier        : error return code
c		      ier=0	returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c---------------------------------------------------------------------
      implicit none
      integer iffld,nterms
      integer ier,ipp,lpp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused
      real *8 rscale,center(3),ztarg(3)
      real *8, allocatable :: w(:)
      complex *16 pot,fld(3)
      complex *16 locexp(0:nterms,-nterms:nterms)
c
c ... Assigning work spaces for various temporary arrays:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+3
      ippd  = ipp+lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr= nterms+3
c
      ifrder=ifr+lfr
      lfrder=nterms+3
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      call l3dtaeval0(rscale,center,locexp,nterms,ztarg,
     1	     pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifr),
     2       w(ifrder))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaeval0(rscale,center,locexp,nterms,
     1		ztarg,pot,iffld,fld,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dtaeval for comments.
c     (pp and ppd are storage arrays for Ynm and Ynm')
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,iffld,i,l,m,n
      real *8 rscale,center(3),ztarg(3),zdiff(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      real *8 fruse,fr(0:nterms+1),frder(0:nterms+1)
      real *8 done,r,theta,phi,d,ctheta,stheta,cphi,sphi
      real *8 phix,phiy,phiz,rs,rx,ry,rz
      real *8 thetax,thetay,thetaz
      complex *16 pot,fld(3),ephi1,ephi1inv
      complex *16 locexp(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
c
      complex *16 eye,ur,utheta,uphi
      complex *16 ztmp,z
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum
      complex *16 ux,uy,uz
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
c     Convert to spherical coordinates
c
      call cart2polarl(zdiff,r,theta,phi)
      d = rscale*r
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute e^{eye*m*phi} array.
c
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i) = i*fr(i-1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     We compensate for this omission by using one lower power in the 
c     r variable - see fruse below.
c     For the n=0 mode, it is not relevant. 
c     
c
      if (iffld.eq.1) then
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
         call ylgndr2s(nterms,ctheta,pp,ppd)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               pp(l,m) = pp(l,m)*rs
               ppd(l,m) = ppd(l,m)*rs
            enddo
         enddo
      else
         call ylgndr(nterms,ctheta,pp)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               pp(l,m) = pp(l,m)*rs
            enddo
         enddo
      endif
c
c
      pot=locexp(0,0)*fr(0)
      if (iffld.eq.1) then
         ur = 0.0d0
         utheta = 0.0d0
         uphi = 0.0d0
      endif
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms
            pot=pot+locexp(n,0)*fr(n)*pp(n,0)
	    ur = ur + frder(n)*pp(n,0)*locexp(n,0)
	    fruse = fr(n-1)*rscale
	    utheta = utheta -locexp(n,0)*fruse*ppd(n,0)*stheta
	    do m=1,n
	       ztmp1=fr(n)*pp(n,m)*stheta
	       ztmp2 = locexp(n,m)*ephi(m) 
	       ztmp3 = locexp(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + frder(n)*pp(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fruse*ppd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fruse*pp(n,m)*ztmpsum
            enddo
         enddo
ccc	 call prin2(' ur is *',ur,2)
ccc	 call prin2(' utheta is *',utheta,2)
ccc	 call prin2(' uphi is *',uphi,2)
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms
	    pot=pot+locexp(n,0)*fr(n)*pp(n,0)
	    do m=1,n
	       ztmp1=fr(n)*pp(n,m)
	       ztmp2 = locexp(n,m)*ephi(m)+locexp(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta(ier,rscale,sources,charge,ns,center,
     1		           nterms,locexp)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  deprecated but left in calling sequence for
c		  backward compatibility.
c
c     locexp    : coeffs for the j-expansion
c
c----------------------------------------------------------------------
      implicit none
      integer ns,l,m,i,ier,nterms
      real *8 rscale,sources(3,ns),center(3),rs
      complex *16 locexp(0:nterms,-nterms:nterms), charge(ns)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call l3dformta1(ier,rscale,sources(1,i),charge(i),
     1		center,nterms,locexp)
      enddo
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            locexp(l,m) = locexp(l,m)*rs
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
c
c**********************************************************************
      subroutine l3dformta1(ier,rscale,source,charge,center,
     &		nterms,locexp)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single charge located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c INPUT:
c
c     rscale    : scaling parameter
c     source    : coordinates of the source
c     charge    : coordinates of the source
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c---------------------------------------------------------------------
c OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   deprecated but left in calling sequence for
c		   backward compatibility.
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,iephi,lephi,ifr,lfr,lused
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      complex *16 locexp(0:nterms,-nterms:nterms), charge
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      lused=ifr+lfr
      allocate(w(lused))
c
      call l3dformta0(rscale,source,charge,center,
     &		nterms,locexp,w(ipp),w(iephi),w(ifr))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0(rscale,source,charge,
     &		center,nterms,locexp,pp,ephi,fr)
c**********************************************************************
c
c     See l3dformta/l3dformta1 for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,n,m,i
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 done,r,theta,phi
      real *8 ctheta,stheta,cphi,sphi,d
      complex *16 fr(0:nterms+1)
      complex *16 locexp(0:nterms,-nterms:nterms), charge
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 ztmp,z
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
c
c     get the Ynm
c
      call ylgndr(nterms,ctheta,pp)
c
c     compute radial functions and scale them by charge strength.
c
      do n = 0, nterms
         fr(n) = fr(n)*charge
      enddo
c
c     Compute contributions to locexp
c
      locexp(0,0)=locexp(0,0) + fr(0)
      do n=1,nterms
         locexp(n,0)=locexp(n,0) + pp(n,0)*fr(n)
         do m=1,n
            ztmp=pp(n,m)*fr(n)
	    locexp(n,m)=locexp(n,m) + ztmp*ephi(-m)
	    locexp(n,-m)=locexp(n,-m) + ztmp*ephi(m)
         enddo
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotfld3dall(iffld,sources,charge,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of charges at 
c     SOURCE(3,ns). 
c     
c       pot =  sum_{i=1,..,ns}  charge(i)/| target - sources(*,i)|
c	fld =  -grad(pot)
c
c     It calls a subroutine for each source.
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     sources(3,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (real *8)        : calculated potential
c     fld   (real *8)        : calculated gradient
c
c---------------------------------------------------------------------
      implicit none
      integer i,ns,iffld
      real *8 sources(3,ns),target(3)
      complex *16 pot,fld(3),potloc,fldloc(3)
      complex *16 eye
      complex *16 charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotfld3d(iffld,sources(1,i),charge(i),target,
     1        potloc,fldloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_targ(iffld,sources,charge,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of charges at 
c     SOURCE(3,ns). 
c     
c       pot =  sum_{i=1,..,ns}  charge(i)/| target - sources(*,i)|
c	fld =  -grad(pot)
c
c     It uses a single loop over all sources.
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     sources(3,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (real *8)        : calculated potential
c     fld   (real *8)        : calculated gradient
c
c---------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      real *8 sources(3,ns),target(3)
      real *8 d,dd,xdiff,ydiff,zdiff,dinv,dinv2,dinv3
      complex *16 pot,fld(3),potloc,fldloc(3)
      complex *16 eye
      complex *16 charge(ns),cd
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if( iffld .eq. 0 ) then
      do i = 1,ns
c
        xdiff=target(1)-sources(1,i)
        ydiff=target(2)-sources(2,i)
        zdiff=target(3)-sources(3,i)
        dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
        d=sqrt(dd)
c
        dinv=1.0d0/d
        pot=pot+charge(i)*dinv
c
      enddo
      endif
c
      if( iffld .eq. 1 ) then
      do i = 1,ns
c
        xdiff=target(1)-sources(1,i)
        ydiff=target(2)-sources(2,i)
        zdiff=target(3)-sources(3,i)
        dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
        d=sqrt(dd)
c
        dinv=1.0d0/d
        pot=pot+charge(i)*dinv
c
ccc        if (iffld.eq.1) then
        dinv2=dinv*dinv
        dinv3=dinv*dinv2
        fld(1)=fld(1)+charge(i)*xdiff*dinv3
        fld(2)=fld(2)+charge(i)*ydiff*dinv3
        fld(3)=fld(3)+charge(i)*zdiff*dinv3
ccc        endif
c
      enddo
      endif
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3d(iffld,source,charge,target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a charge at 
c     SOURCE. 
c     
c              	pot = charge/ |target-source|
c		fld = -grad(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld     : flag for computing gradient
c	                 	iffld = 0 -> dont compute 
c		                iffld = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     fld       : calculated gradient
c
c---------------------------------------------------------------------
      implicit none
      integer iffld
      real *8 source(3),target(3)
      real *8 xdiff,ydiff,zdiff,dd,d,dinv,dinv2,dinv3
      complex *16 pot,fld(3)
      complex *16 h0,h1,cd,eye,z,ewavek
      complex *16 charge
c
      data eye/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      dinv=1.0d0/d
      pot=charge*dinv
c
      if (iffld.eq.1) then
         dinv2=dinv*dinv
         dinv3=dinv*dinv2
         fld(1)=charge*xdiff*dinv3
         fld(2)=charge*ydiff*dinv3
         fld(3)=charge*zdiff*dinv3
      endif
      return
      end
c
c**********************************************************************
      subroutine l3dadd(mpole,mpole2,nterms)
c**********************************************************************
c
c     add mpole to mpole2
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,i,j
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
c
      do i = 0,nterms
         do j = -i,i
	    mpole2(i,j) = mpole2(i,j)+mpole(i,j)
	 enddo
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dadd_trunc(mpole,mpole2,nterms,ldc)
c**********************************************************************
c
c     add mpole to mpole2, assuming size of mpole is smaller than
c     size of mpole2 (nterms < ldc).
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,ldc,i,j
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:ldc,-ldc:ldc)
c
      do i = 0,nterms
         do j = -i,i
	    mpole2(i,j) = mpole2(i,j)+mpole(i,j)
	 enddo
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine cart2polarl(zat,r,theta,phi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c	zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c	r     :  |zat|
c	theta : angle subtended with respect to z-axis
c	phi   : angle of (zat(1),zat(2)) subtended with 
c               respect to x-axis
c
c-----------------------------------------------------------------------
      implicit none
      real *8 zat(3),r,proj,theta,phi
      complex *16 ephi,eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2+zat(3)**2)
      proj = sqrt(zat(1)**2+zat(2)**2)
c
      theta = datan2(proj,zat(3))
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      phi = 0
      else
      phi = datan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
c
c
c**********************************************************************
        subroutine l3drhpolar(x,y,z,r,ctheta,ephi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       x,y,z   : Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r      : sqrt(x*x+y*y+z*z)
c       ctheta : cos(theta)
c       ephi   : exp(I*phi)  (complex *16_
c
c       where
c
c       theta is angle subtended with respect to z-axis
c       phi   is angle of (x,y) subtended with 
c               respect to x-axis
c-----------------------------------------------------------------------
        implicit none
        real *8 x,y,z,r,ctheta,proj
        complex *16 ephi,ima
        data ima/(0.0d0,1.0d0)/
c
        proj = sqrt(x*x+y*y)
        r = sqrt(x*x+y*y+z*z)
c
        if( abs(r) .gt. 0 ) then
        ctheta = z/r
        else
        ctheta = 0.0d0
        endif
c
        if( abs(proj) .gt. 0 ) then
        ephi = cmplx(x,y)/proj
        else
        ephi = 0.0d0
        endif
c
        return
        end
c
c
c
c
c
C***********************************************************************
      subroutine l3dformmp_dp(ier,rscale,sources,dipstr,dipvec,ns,
     1                  center,nterms,mpole)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS 
c     dipole sources C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipstr(ns)      : source strengths
C     dipvec(3,ns)    : dipole vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		         ier=0  returned successfully
c		         deprecated but left in calling sequence for
c		         backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c                  
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused
      real *8 center(3),sources(3,ns)
      real *8 dipvec(3,ns)
      real *8 rscale,proj,rs
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,dipstr(ns)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1, ns
         call l3dformmp1_dp(ier,rscale,sources(1,i),dipstr(i),
     1        dipvec(1,i),center,nterms,mpole)
      enddo
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            mpole(l,m) = mpole(l,m)*rs
         enddo
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_dp(ier,rscale,source,dipstr,dipvec,
     1		center,nterms,mpole)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a dipole located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     dipstr  : complex dipole strength
c     dipvec  : dipole direction vector
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,ippd,iephi,lephi,ifrder,lfrder,ifr,lused,jer
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 dipvec(3)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 dipstr
c
c     compute workspace requirements
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifrder=iephi+lephi
      lfrder=2*(nterms+3)
c
      ifr=ifrder+lfrder
      lused=ifr + lfrder
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp0_dp(jer,rscale,source,dipstr,dipvec,
     1		center,nterms,mpole,w(ipp),w(ippd),w(iephi),
     2          w(ifr),w(ifrder))
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_dp(ier,rscale,source,dipstr,dipvec,
     1		center,nterms,mpole,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dformmp1_dp for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,i,nterms,n,m
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      real *8 r,theta,phi,d
      real *8 ctheta,stheta,cphi,sphi
      real *8 phix,phiy,phiz
      real *8 rx,ry,rz,thetax,thetay,thetaz
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 dipstr
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fr(0:nterms+1),ztmp,frder(0:nterms+1),z
      complex *16 fruse,ux,uy,uz,ur,utheta,uphi,zzz
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = r
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      d = d*rscale
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i) = i*fr(i-1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     For the n=0 mode, it is not relevant. For n>0 modes,
c     the variable fruse is set to fr(n)/r:
c
c     
c
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
c
c     get the associated Legendre functions:
c
      call ylgndr2s(nterms,ctheta,pp,ppd)
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r = 
c         \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)/ (|T|^(n+1))
c
c     so contribution is |S|^n times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      ur = pp(0,0)*frder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      ux = ur*rx + utheta*thetax + uphi*phix
      uy = ur*ry + utheta*thetay + uphi*phiy
      uz = ur*rz + utheta*thetaz + uphi*phiz
      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      mpole(0,0)= mpole(0,0) + zzz*dipstr
      do n=1,nterms
         fruse = fr(n-1)*rscale
         ur = pp(n,0)*frder(n)
         utheta = -fruse*ppd(n,0)*stheta
         uphi = 0.0d0
         ux = ur*rx + utheta*thetax + uphi*phix
         uy = ur*ry + utheta*thetay + uphi*phiy
         uz = ur*rz + utheta*thetaz + uphi*phiz
         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         mpole(n,0)= mpole(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fruse*ppd(n,m)
            uphi = -eye*m*ephi(-m)*fruse*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            mpole(n,m)= mpole(n,m) + zzz*dipstr
c
            ur = frder(n)*pp(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fruse*ppd(n,m)
            uphi = eye*m*ephi(m)*fruse*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            mpole(n,-m)= mpole(n,-m) + zzz*dipstr
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_dp(iffld,sources,dipstr,dipvec,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of dipoles at 
c     SOURCE(3,ns). 
c     
c     The potential due to a single dipole is 
c
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	 fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                   iffld = 0 -> dont compute 
c		           iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      real *8 sources(3,ns),target(3)
      real *8 dipvec(3,ns)
      complex *16 pot,fld(3),potloc,fldloc(3)
      complex *16 eye
      complex *16 dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotfld3d_dp(iffld,sources(1,i),dipstr(i),dipvec(1,i),
     1        target,potloc,fldloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_dp_targ(iffld,sources,dipstr,dipvec,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of dipoles at 
c     SOURCE(3,ns). 
c     
c     The potential due to a single dipole is 
c
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	 fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                   iffld = 0 -> dont compute 
c		           iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      real *8 sources(3,ns),target(3)
      real *8 dipvec(3,ns)
      complex *16 pot,fld(3),potloc,fldloc(3)
      complex *16 eye
      complex *16 dipstr(ns)
      real *8 xdiff,ydiff,zdiff,dd,d,dinv,dinv2,
     $   dinv3,ddd,dotprod,dinv5,rtttt
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if( iffld .eq. 0 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)
      pot=pot+ dipstr(i)*(dotprod*dinv3)
c
      enddo
      endif
c
      if( iffld .eq. 1 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)
      pot=pot+ dipstr(i)*(dotprod*dinv3)
c
ccc      if (iffld.eq.1) then
         dinv5 = dinv3*dinv2
         rtttt = 3.0d0*dotprod*dinv5
         fld(1)=fld(1)+dipstr(i)*(rtttt*xdiff-dinv3*dipvec(1,i))
         fld(2)=fld(2)+dipstr(i)*(rtttt*ydiff-dinv3*dipvec(2,i))
         fld(3)=fld(3)+dipstr(i)*(rtttt*zdiff-dinv3*dipvec(3,i))
ccc      endif 
c
      enddo
      endif
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3d_dp(iffld,source,dipstr,dipvec,target,
     1                        pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     source(3)    : location of the source 
c     dipstr(ns)   : dipole strength
c     dipvec(3,ns) : dipole direction
c     target(3)    : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot          : calculated potential
c     fld          : calculated -gradient
c
c----------------------------------------------------------------------
      implicit none
      integer iffld
      real *8 source(3),target(3)
      real *8 dipvec(3)
      real *8 xdiff,ydiff,zdiff,dd,d,dinv,dinv2,
     $   dinv3,ddd,dotprod,dinv5,rtttt
      complex *16 pot,fld(3)
      complex *16 cd,eye,ztttt,cd2
      complex *16 dipstr,z1,z2,z3
c
      data eye/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
c ... Calculate the potential and field in the regular case:
c
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)+zdiff*dipvec(3)
      pot= dipstr*(dotprod*dinv3)
      if (iffld.eq.1) then
         dinv5 = dinv3*dinv2
         rtttt = 3.0d0*dotprod*dinv5
         fld(1)=dipstr*(rtttt*xdiff-dinv3*dipvec(1))
         fld(2)=dipstr*(rtttt*ydiff-dinv3*dipvec(2))
         fld(3)=dipstr*(rtttt*zdiff-dinv3*dipvec(3))
      endif 
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta_dp(ier,rscale,sources,dipstr,dipvec,ns,
     1		           center,nterms,locexp)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  deprecated but left in calling sequence for
c		  backward compatibility.
c
c     locexp    : coeffs for the j-expansion
c
c
c----------------------------------------------------------------------
      implicit none
      integer ier,ns,nterms,l,m,i
      real *8 rscale,sources(3,ns),center(3),rs
      real *8 dipvec(3,ns)
      complex *16 locexp(0:nterms,-nterms:nterms), dipstr(ns)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call l3dformta1_dp(ier,rscale,sources(1,i),dipstr(i),
     1		dipvec(1,i),center,nterms,locexp)
      enddo
c
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            locexp(l,m) = locexp(l,m)*rs
         enddo
      enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta1_dp(ier,rscale,source,dipstr,dipvec,
     &		center,nterms,locexp)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single dipole located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   deprecated but left in calling sequence for
c		   backward compatibility.
c
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 dipvec(3)
      complex *16 locexp(0:nterms,-nterms:nterms), dipstr
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      ippd = ipp+lpp
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      ifrder=ifr+lfr
      lfrder=2*(nterms+3)
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      call l3dformta0_dp(rscale,source,dipstr,dipvec,
     &   center,nterms,locexp,w(ipp),w(ippd),w(iephi),w(ifr),w(ifrder))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0_dp(rscale,source,dipstr,dipvec,
     &		center,nterms,locexp,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dformta_dp/l3dformta1_dp for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,i,n,m
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      real *8 r,theta,phi,ctheta,stheta,cphi,sphi,done
      real *8 d,phix,phiy,phiz,rx,ry,rz,thetax,thetay,thetaz
      complex *16 locexp(0:nterms,-nterms:nterms), dipstr
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fr(0:nterms+1),ztmp,frder(0:nterms+1),z
      complex *16 ux,uy,uz,ur,utheta,uphi,zzz
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
      fr(nterms+1)=fr(nterms)*d
      do i=0,nterms
         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
        rx = stheta*cphi
        thetax = ctheta*cphi/r
        phix = -sphi/r
        ry = stheta*sphi
        thetay = ctheta*sphi/r
        phiy = cphi/r
        rz = ctheta
        thetaz = -stheta/r
        phiz = 0.0d0
c
c
c     get the associated Legendre functions:
c
      call ylgndr2s(nterms,ctheta,pp,ppd)
c
c     Compute contribution to local coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r = 
c         \sum_n 1/(2n+1) \sum_m  |T|^n Ylm(T) Ylm*(S) / (|S|^{n+1})
c
c     so contribution is |S|^{n+1} times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      ur = pp(0,0)*frder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      ux = ur*rx + utheta*thetax + uphi*phix
      uy = ur*ry + utheta*thetay + uphi*phiy
      uz = ur*rz + utheta*thetaz + uphi*phiz
      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      locexp(0,0)= locexp(0,0) + zzz*dipstr
      do n=1,nterms
         ur = pp(n,0)*frder(n)
         utheta = -fr(n)*ppd(n,0)*stheta
         uphi = 0.0d0
         ux = ur*rx + utheta*thetax + uphi*phix
         uy = ur*ry + utheta*thetay + uphi*phiy
         uz = ur*rz + utheta*thetaz + uphi*phiz
         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         locexp(n,0)= locexp(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fr(n)*ppd(n,m)
            uphi = -eye*m*ephi(-m)*fr(n)*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            locexp(n,m)= locexp(n,m) + zzz*dipstr
c
            ur = frder(n)*pp(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fr(n)*ppd(n,m)
            uphi = eye*m*ephi(m)*fr(n)*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            locexp(n,-m)= locexp(n,-m) + zzz*dipstr
         enddo
      enddo
c
      return
      end
c
c
c**********************************************************************
      subroutine lpotfld3dall_sdp_targ(iffld,sources,
     $     charge,dipstr,dipvec,ns,target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection 
c     of charges and dipoles at SOURCE(3,ns). 
c     
c     The potential due to a single charge is 
c
c        pot = charge/r 
c
c     and the potential due to a single dipole is 
c
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	 fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                   iffld = 0 -> dont compute 
c		           iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     charge(ns)    : charge strength
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      real *8 sources(3,ns),target(3)
      real *8 dipvec(3,ns)
      complex *16 pot,fld(3),potloc,fldloc(3)
      complex *16 eye
      complex *16 charge(ns),dipstr(ns)
      real *8 xdiff,ydiff,zdiff,dd,d,dinv,dinv2,
     $   dinv3,ddd,dotprod,dinv5,rtttt
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if( iffld .eq. 0 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)

      pot=pot+charge(i)*dinv
      pot=pot+dipstr(i)*(dotprod*dinv3)
c
      enddo
      endif
c
      if( iffld .eq. 1 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)
c
      pot=pot+charge(i)*dinv
      pot=pot+dipstr(i)*(dotprod*dinv3)
c
ccc      if (iffld.eq.1) then
         dinv5 = dinv3*dinv2
         rtttt = 3.0d0*dotprod*dinv5
         fld(1)=fld(1)+charge(i)*xdiff*dinv3
         fld(2)=fld(2)+charge(i)*ydiff*dinv3
         fld(3)=fld(3)+charge(i)*zdiff*dinv3
         fld(1)=fld(1)+dipstr(i)*(rtttt*xdiff-dinv3*dipvec(1,i))
         fld(2)=fld(2)+dipstr(i)*(rtttt*ydiff-dinv3*dipvec(2,i))
         fld(3)=fld(3)+dipstr(i)*(rtttt*zdiff-dinv3*dipvec(3,i))
ccc      endif 
c
      enddo
      endif
c
      return
      end
c
c
c
c
