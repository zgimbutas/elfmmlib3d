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
c      This file contains the accelerated subroutines for 
c      forming and evaluating multipole expansions.
c
c
c**********************************************************************
      subroutine l3dmpevalall_trunc(rscale,center,mpole,nterms,nterms1,
     $     ztarg,nt,ifpot,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to a TRUNCATED outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld = -gradient(pot) if iffld = 1.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :   scaling parameter (see formmp1l3d)
c     center :   expansion center
c     mpole  :   multipole expansion in 2d matrix format
c     nterms :   order of the multipole expansion
c     nterms1 :   order of truncated expansion to be used
c     ztarg  :   target location
c     nt     :   number of targets
c     ifpot  :   flag controlling evaluation of potential
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.        
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.        
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     fld    :    gradient at ztarg (if requested)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lused,i
      real *8 rscale,center(3),ztarg(3,nt)
      real *8 wlege(0:nlege,0:nlege)
      real *8, allocatable :: w(:)
      complex *16 pot(nt),fld(3,nt)
      complex *16 pot0,fld0(3)
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
      do i=1,nt
      call l3dmpeval_trunc0(rscale,center,mpole,nterms,nterms1,
     $   ztarg(1,i),pot0,iffld,fld0,w(ipp),w(ippd),
     $     w(iephi),w(ifr),w(ifrder),wlege,nlege)
      if( ifpot .eq. 1 ) pot(i)=pot(i)+pot0
      if( iffld .eq. 1 ) then
        fld(1,i)=fld(1,i)+fld0(1)
        fld(2,i)=fld(2,i)+fld0(2)
        fld(3,i)=fld(3,i)+fld0(3)
      endif
      enddo
c
ccc      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dmpeval_trunc(rscale,center,mpole,nterms,nterms1,
     $     ztarg,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
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
c     nterms1 :   order of truncated expansion to be used
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lused,i
      real *8 rscale,center(3),ztarg(3)
      real *8 wlege(0:nlege,0:nlege)
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
      call l3dmpeval_trunc0(rscale,center,mpole,nterms,nterms1,ztarg,
     1	   pot,iffld,fld,w(ipp),w(ippd),
     $     w(iephi),w(ifr),w(ifrder),wlege,nlege)
ccc      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dmpeval_trunc0(rscale,center,mpole,nterms,nterms1,
     1		ztarg,pot,iffld,fld,ynm,ynmd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dmpeval for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,iffld,ifpot,nlege
      integer n,m,i,l
      real *8 rscale,center(3),ztarg(3),zdiff(3)
      real *8 wlege(0:nlege,0:nlege)
      real *8 ynm(0:nterms1,0:nterms1)
      real *8 ynmd(0:nterms1,0:nterms1)
      real *8 fr(0:nterms+1)
      real *8 frder(0:nterms+1)
      real *8 done,r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      complex *16 pot,fld(3),ephi1,ur,utheta,uphi,ux,uy,uz
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
c
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
      ctheta = dcos(theta)
      d = 1.0d0/r
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      d = 1/r
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
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
         ephi(-i)=conjg(ephi(i))
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
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
ccc         call ylgndr2s(nterms1,ctheta,ynm,ynmd)
c         call ylgndr2sfw(nterms1,ctheta,ynm,ynmd,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               ynm(l,m) = ynm(l,m)*rs
c               ynmd(l,m) = ynmd(l,m)*rs
c            enddo
c         enddo
         call ylgndru2sfw(nterms1,ctheta,ynm,ynmd,wlege,nlege)
      else        
ccc         call ylgndr(nterms1,ctheta,ynm)
c         call ylgndrfw(nterms1,ctheta,ynm,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               ynm(l,m) = ynm(l,m)*rs
c            enddo
c         enddo
         call ylgndrufw(nterms1,ctheta,ynm,wlege,nlege)
      endif
c
c
c
c     initialize computed values and 
c     scale derivatives of Hankel functions so that they are
c     derivatives with respect to r.
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
         do n=1,nterms1
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
         do n=1,nterms1
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
c**********************************************************************
      subroutine l3dtaevalall_trunc(rscale,center,locexp,nterms,nterms1,
     1		ztarg,nt,ifpot,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates a TRUNCATED local expansion centered 
c     at CENTER at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) r^n Y_nm(theta,phi)
c             n   m
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale     : scaling parameter used in forming expansion
c                                   (see l3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the local expansion
c     nterms     : order of the local expansion
c     nterms1    : order of the truncated expansion to be used
c     ztarg      : vector of targets 
c     nt         : number of targets
c     ifpot      : flag for potential computation
c		                    ifpot=0  - pot is not computed
c		                    ifpot=1  - pot is computed
c     iffld      : flag for gradient computation
c		                    iffld=0  - gradient is not computed
c		                    iffld=1  - gradient is computed
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot        : potential aat ztarg (if requested)
c     fld(3)     : gradient at ztarg (if requested)
c     ier        : error return code
c		      ier=0	returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused,i
      real *8 rscale,center(3),ztarg(3,nt)
      real *8 wlege(0:nlege,0:nlege)
      real *8, allocatable :: w(:)
      complex *16 pot(nt),fld(3,nt)
      complex *16 pot0,fld0(3)
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
      do i=1,nt
      call l3dtaeval_trunc0(rscale,center,locexp,nterms,nterms1,
     $   ztarg(1,i),pot0,iffld,fld0,w(ipp),w(ippd),w(iephi),w(ifr),
     2   w(ifrder),wlege,nlege)
      if( ifpot .eq. 1 ) pot(i)=pot(i)+pot0
      if( iffld .eq. 1 ) then
        fld(1,i)=fld(1,i)+fld0(1)
        fld(2,i)=fld(2,i)+fld0(2)
        fld(3,i)=fld(3,i)+fld0(3)
      endif
      enddo
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaeval_trunc(rscale,center,locexp,nterms,nterms1,
     1		ztarg,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates a TRUNCATED local expansion centered 
c     at CENTER at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) r^n Y_nm(theta,phi)
c             n   m
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale     : scaling parameter used in forming expansion
c                                   (see l3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the local expansion
c     nterms     : order of the local expansion
c     nterms     : order of the truncated expansion to be used
c     ztarg      : target vector
c     iffld      : flag for gradient computation
c		                    iffld=0  - gradient is not computed
c		                    iffld=1  - gradient is computed
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot        : potential at ztarg(3)
c     fld(3)     : gradient at ztarg (if requested)
c     ier        : error return code
c		      ier=0	returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,iffld,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused,i
      real *8 rscale,center(3),ztarg(3)
      real *8 wlege(0:nlege,0:nlege)
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
      call l3dtaeval_trunc0(rscale,center,locexp,nterms,nterms1,ztarg,
     1	     pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifr),
     2       w(ifrder),wlege,nlege)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaeval_trunc0(rscale,center,locexp,nterms,nterms1,
     1		ztarg,pot,iffld,fld,pp,ppd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dtaeval for comments.
c     (pp and ppd are storage arrays for Ynm and Ynm')
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,iffld,nlege
      integer i,l,n,m
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale,center(3),ztarg(3),zdiff(3)
      real *8 pp(0:nterms1,0:nterms1)
      real *8 ppd(0:nterms1,0:nterms1)
      real *8 fruse,fr(0:nterms+1),frder(0:nterms+1)
      real *8 done,r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
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
      ctheta = dcos(theta)
      d = rscale*r
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      d = rscale*r
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
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
         ephi(-i)=conjg(ephi(i))
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
c     For the n=0 mode, it is not relevant. 
c
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
c         call ylgndr2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               pp(l,m) = pp(l,m)*rs
c               ppd(l,m) = ppd(l,m)*rs
c            enddo
c         enddo
         call ylgndru2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
      else
c         call ylgndrfw(nterms1,ctheta,pp,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               pp(l,m) = pp(l,m)*rs
c            enddo
c         enddo
         call ylgndrufw(nterms1,ctheta,pp,wlege,nlege)
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
         do n=1,nterms1
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
         do n=1,nterms1
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
C***********************************************************************
      subroutine l3dformmp_trunc(ier,rscale,sources,charge,ns,center,
     1                  nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     nterms1         : order of truncated expansion to be generated
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		         ier=0  returned successfully
c    		         deprecated but left in calling sequence for
c		         backward compatibility.
c    
c     mpole           : coeffs of the multipole expansion
c                  
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,nlege
      integer ns,i,l,m,ier1
      real *8 center(3),sources(3,ns)
      real *8 rscale,rs
      real *8 wlege(0:nlege,0:nlege)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,charge(ns)
      integer ipp,lpp,iephi,lephi,ifr,lfr,lused
      real *8, allocatable :: w(:)
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
c
c ... Assign work spaces:
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
      do i = 1, ns
c         call l3dformmp_trunc1
c     $   (ier1,rscale,sources(1,i),charge(i),center,
c     1        nterms,nterms1,mpole,wlege,nlege)
        call l3dformmp_trunc0(rscale,sources(1,i),charge(i),center,
     $   nterms,nterms1,
     1   mpole,w(ipp),w(iephi),w(ifr),wlege,nlege)
      enddo
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            mpole(l,m) = mpole(l,m)*rs
c         enddo
c      enddo
c
      return
      end
C
C***********************************************************************
      subroutine l3dformmp_add_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs TRUNCATED multipole expansion about CENTER due to 
C     NS sources located at SOURCES(3,*)  and INCREMENTS mpole
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
C     nterms1         : order of truncated multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		           ier=0  returned successfully
c		           deprecated but left in calling sequence for
c		           backward compatibility.
c    
c     mpole           : incremented multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer ns,i,l,m,ier,ier1
      real *8 center(3),sources(3,ns)
      real *8 rscale
      real *8 wlege(0:nlege,0:nlege)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,charge(ns)
      complex *16, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
        allocate( mptemp(0:nterms,-nterms:nterms) )
C
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo
c
        call l3dformmp_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mptemp,wlege,nlege)
c
        do l = 0,nterms
          do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
          enddo
        enddo
c
      return
      end
c
C
c**********************************************************************
      subroutine l3dformmp_trunc1(ier,rscale,source,charge,center,
     1		nterms,nterms1,mpole,wlege,nlege)
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
c     nterms  : order of the multipole expansion
c     nterms1 : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c                            
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,lused,i
      real *8 rscale,source(3),center(3)
      real *8 wlege(0:nlege,0:nlege)
      real *8, allocatable :: w(:)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 charge
c
c ... Assign work spaces:
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
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp_trunc0(rscale,source,charge,center,
     $   nterms,nterms1,
     1   mpole,w(ipp),w(iephi),w(ifr),wlege,nlege)
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp_trunc0(rscale,source,charge,center,
     1		nterms,nterms1,mpole,pp,ephi,fr,wlege,nlege)
c**********************************************************************
c
c     See l3dformmp1 for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer n,m,i
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 wlege(0:nlege,0:nlege)
      real *8 pp(0:nterms1,0:nterms1)
      real *8 done,r,theta,phi,dtmp
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 charge
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16  fr(0:nterms+1)
      complex *16  ztmp,z
c
c
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      d = r
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      d = r
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
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
         ephi(-i)=conjg(ephi(i))
      enddo
c
c     get the associated Legendre functions:
c
ccc      call ylgndrfw(nterms1,ctheta,pp,wlege,nlege)
      call ylgndrufw(nterms1,ctheta,pp,wlege,nlege)
ccc      call prinf(' after ylgndr with nterms = *',nterms,1)
ccc      call prinm2(pp,nterms)
c
c     multiply all fr's by charge strength.
c
      do n = 0,nterms1
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
      do n=1,nterms1
         mpole(n,0)= mpole(n,0) + pp(n,0)*fr(n)
         do m=1,n
cc            ztmp=pp(n,m)*fr(n)
cc            mpole(n, m)= mpole(n, m) + ztmp*dconjg(ephi(m))
cc            mpole(n,-m)= mpole(n,-m) + ztmp*dconjg(ephi(-m))
            ztmp=pp(n,m)*fr(n)
            mpole(n, m)= mpole(n, m) + ztmp*ephi(-m)
            mpole(n,-m)= mpole(n,-m) + ztmp*ephi(+m)
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
      subroutine l3dformta_trunc(ier,rscale,sources,charge,ns,center,
     1		        nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
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
c     nterms    : order of the local expansion
c     nterms1    : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : coeffs for the local expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer ns,i,l,m
      real *8 rscale,rs,sources(3,ns),center(3)
      real *8 wlege(0:nlege,0:nlege)
      complex *16 locexp(0:nterms,-nterms:nterms), charge(ns)
      complex *16 eye
      integer ipp,lpp,iephi,lephi,ifr,lfr,lused
      real *8, allocatable :: w(:)
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
      do i = 1,ns
c         call l3dformta_trunc1(ier,rscale,sources(1,i),charge(i),
c     1		center,nterms,nterms1,locexp,wlege,nlege)
         call l3dformta_trunc0(rscale,sources(1,i),charge(i),center,
     &      nterms,nterms1,locexp,w(ipp),w(iephi),w(ifr),
     $      wlege,nlege)
      enddo
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            locexp(l,m) = locexp(l,m)*rs
c         enddo
c      enddo
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta_add_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below. INCREMENT
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : incremented local expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer ns,l,m
      real *8 rscale,sources(3,ns),center(3)
      real *8 wlege(0:nlege,0:nlege)
      complex *16 locexp(0:nterms,-nterms:nterms), charge(ns)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo
c
      call l3dformta_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m) + mptemp(l,m)
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
      subroutine l3dformta_trunc1(ier,rscale,source,charge,center,
     &		nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single charge located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c     source    : coordinates of the source
c     charge    : coordinates of the source
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   deprecated but left in calling sequence for
c		   backward compatibility.
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,lused
      real *8 rscale,source(3),center(3)
      real *8 wlege(0:nlege,0:nlege)
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
      call l3dformta_trunc0(rscale,source,charge,center,
     &   nterms,nterms1,locexp,w(ipp),w(iephi),w(ifr),
     $   wlege,nlege)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta_trunc0(rscale,source,charge,
     &		center,nterms,nterms1,locexp,pp,ephi,fr,wlege,nlege)
c**********************************************************************
c
c     See l3dformta/l3dformta1 for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer i,n,m
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 pp(0:nterms1,0:nterms1)
      real *8 wlege(0:nlege,0:nlege)
      real *8 done,r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
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
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=conjg(ephi1)
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=conjg(ephi(i))
      enddo
c
c     get the Ynm
c
      call ylgndrufw(nterms1,ctheta,pp,wlege,nlege)
c
c     compute radial functions and scale them by charge strength.
c
      do n = 0, nterms1
         fr(n) = fr(n)*charge
      enddo
c
c     Compute contributions to locexp
c
      locexp(0,0)=locexp(0,0) + fr(0)
      do n=1,nterms1
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
c
C***********************************************************************
      subroutine l3dformmp_dp_trunc(ier,rscale,sources,dipstr,dipvec,ns,
     1                  center,nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs TRUNCATED multipole expansion about CENTER due to NS 
c     dipole sources located at SOURCES(3,*).
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
C     nterms1         : order of truncated multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		          ier=0  returned successfully
c		          deprecated but left in calling sequence for
c		          backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ns,i,l,m, ier,nlege
      real *8 center(3),sources(3,ns)
      real *8 dipvec(3,ns)
      real *8 rscale,rs
      real *8 wlege(0:nlege,0:nlege)
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
         call l3dformmp1_dp_trunc(ier,rscale,sources(1,i),dipstr(i),
     1        dipvec(1,i),center,nterms,nterms1,mpole,wlege,nlege)
      enddo
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            mpole(l,m) = mpole(l,m)*rs
c         enddo
c      enddo
c
      return
      end
C
C***********************************************************************
      subroutine l3dformmp_dp_add_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     dipole sources located at SOURCES(3,*) and INCREMENTS mpole.
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
C     nterms1         : order of truncated multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		          ier=0  returned successfully
c		          deprecated but left in calling sequence for
c		          backward compatibility.
c
c     mpole           : incremented multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ns,i,l,m, ier,nlege
      real *8 center(3),sources(3,ns)
      real *8 dipvec(3,ns)
      real *8 rscale
      real *8 wlege(0:nlege,0:nlege)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
      complex *16 eye,dipstr(ns)
      data eye/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo

      call l3dformmp_dp_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
         enddo
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_dp_trunc(ier,rscale,source,dipstr,dipvec,
     1		center,nterms,nterms1,mpole,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the truncated multipole expansion 
c     about CENTER due to a dipole located at the point SOURCE.
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
c     nterms  : order of the multipole expansion
c     nterms1 : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
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
      integer nterms,nterms1,ier,jer,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lfrder,lused
      real *8 rscale,source(3),center(3)
      real *8 wlege(0:nlege,0:nlege)
      real *8, allocatable :: w(:)
      real *8 dipvec(3)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 dipstr
c
c ... Assign work spaces:
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
      call l3dformmp0_dp_trunc(jer,rscale,source,dipstr,dipvec,
     1		center,nterms,nterms1,mpole,w(ipp),w(ippd),w(iephi),
     2          w(ifr),w(ifrder),wlege,nlege)
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_dp_trunc(ier,rscale,source,dipstr,dipvec,
     1    center,nterms,nterms1,mpole,pp,ppd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dformmp1_dp for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer i,n,m
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 rfac1, rfac2, rfac3
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      real *8 wlege(0:nlege,0:nlege)
      real *8 done,r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
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
ccc      call ylgndr2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
      call ylgndru2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
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
      rfac1 = dipvec(1)*rx + dipvec(2)*ry + dipvec(3)*rz
      rfac2 = dipvec(1)*thetax + dipvec(2)*thetay + dipvec(3)*thetaz
      rfac3 = dipvec(1)*phix + dipvec(2)*phiy + dipvec(3)*phiz
ccc      ux = ur*rx + utheta*thetax + uphi*phix
ccc      uy = ur*ry + utheta*thetay + uphi*phiy
ccc      uz = ur*rz + utheta*thetaz + uphi*phiz
ccc      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
      mpole(0,0)= mpole(0,0) + zzz*dipstr
c
      do n=1,nterms1
         fruse = fr(n-1)*rscale
         ur = pp(n,0)*frder(n)
         utheta = -fruse*ppd(n,0)*stheta
         uphi = 0.0d0
ccc         ux = ur*rx + utheta*thetax + uphi*phix
ccc         uy = ur*ry + utheta*thetay + uphi*phiy
ccc         uz = ur*rz + utheta*thetaz + uphi*phiz
ccc         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
         mpole(n,0)= mpole(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta
            utheta = -fruse*ppd(n,m)
            uphi = -eye*m*fruse*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            zzz = (rfac1*ur + rfac2*utheta + rfac3*uphi)*ephi(-m)
            mpole(n,m)= mpole(n,m) + zzz*dipstr
c
ccc            ur = frder(n)*pp(n,m)*stheta*ephi(m)
ccc            utheta = -ephi(m)*fruse*ppd(n,m)
ccc            uphi = eye*m*ephi(m)*fruse*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
ccc            zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
            zzz = conjg(zzz) 
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
      subroutine l3dformta_dp_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
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
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : coeffs for the local expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ns,ier,nlege
      integer i,l,m
      real *8 sources(3,ns),center(3),rs,rscale
      real *8 dipvec(3,ns)
      real *8 wlege(0:nlege,0:nlege)
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
         call l3dformta1_dp_trunc(ier,rscale,sources(1,i),dipstr(i),
     1      dipvec(1,i),center,nterms,nterms1,locexp,wlege,nlege)
      enddo
c
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            locexp(l,m) = locexp(l,m)*rs
c         enddo
c      enddo
c
      return
      end
c
c
c**********************************************************************
      subroutine l3dformta_dp_add_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*)
c     and INCREMENTS locexp.
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
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : coeffs for the j-expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,ns,nterms1,ier,nlege
      integer l,m
      real *8 rscale,sources(3,ns),center(3)
      real *8 dipvec(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      complex *16 locexp(0:nterms,-nterms:nterms), dipstr(ns)
      complex *16 eye
      complex *16, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo
c
      call l3dformta_dp_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m)+mptemp(l,m)
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
      subroutine l3dformta1_dp_trunc(ier,rscale,source,dipstr,dipvec,
     &		center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the truncated local expansion about 
c     CENTER due to a single dipole located at SOURCE.
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
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused,i
      real *8 rscale,source(3),center(3)
      real *8 wlege(0:nlege,0:nlege)
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
      call l3dformta0_dp_trunc(rscale,source,dipstr,dipvec,
     &   center,nterms,nterms1,locexp,
     $   w(ipp),w(ippd),w(iephi),w(ifr),w(ifrder),wlege,nlege)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0_dp_trunc(rscale,source,dipstr,dipvec,
     &     center,nterms,nterms1,
     $     locexp,pp,ppd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dformta_dp_trunc/l3dformta1_dp_trunc for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer i,n,m
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 rfac1, rfac2, rfac3
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      real *8 wlege(0:nlege,0:nlege)
      real *8 done,r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
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
ccc      call ylgndr2sfw(nterms,ctheta,pp,ppd,wlege,nlege)
      call ylgndru2sfw(nterms,ctheta,pp,ppd,wlege,nlege)
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
      rfac1 = dipvec(1)*rx + dipvec(2)*ry + dipvec(3)*rz
      rfac2 = dipvec(1)*thetax + dipvec(2)*thetay + dipvec(3)*thetaz
      rfac3 = dipvec(1)*phix + dipvec(2)*phiy + dipvec(3)*phiz
ccc      ux = ur*rx + utheta*thetax + uphi*phix
ccc      uy = ur*ry + utheta*thetay + uphi*phiy
ccc      uz = ur*rz + utheta*thetaz + uphi*phiz
ccc      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
      locexp(0,0)= locexp(0,0) + zzz*dipstr
      do n=1,nterms
         ur = pp(n,0)*frder(n)
         utheta = -fr(n)*ppd(n,0)*stheta
         uphi = 0.0d0
ccc         ux = ur*rx + utheta*thetax + uphi*phix
ccc         uy = ur*ry + utheta*thetay + uphi*phiy
ccc         uz = ur*rz + utheta*thetaz + uphi*phiz
ccc         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
         locexp(n,0)= locexp(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta
            utheta = -fr(n)*ppd(n,m)
            uphi = -eye*m*fr(n)*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            zzz = (rfac1*ur + rfac2*utheta + rfac3*uphi)*ephi(-m)
            locexp(n,m)= locexp(n,m) + zzz*dipstr
c
ccc            ur = frder(n)*pp(n,m)*stheta*ephi(m)
ccc            utheta = -ephi(m)*fr(n)*ppd(n,m)
ccc            uphi = eye*m*ephi(m)*fr(n)*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            zzz = conjg(zzz) 
            locexp(n,-m)= locexp(n,-m) + zzz*dipstr
         enddo
      enddo
c
      return
      end
c
c
C***********************************************************************
c
c
c       Multipole forming routines for real-valued charges and dipoles
c
c
C***********************************************************************
      subroutine l3dformmp_charge_trunc(ier,rscale,sources,charge,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS real-valued
c     charge sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : charge strengths 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      real *8 center(3),sources(3,ns)
      real *8 charge(ns)
      real *8 rscale
      complex *16 mpole(0:nterms,-nterms:nterms)
      integer nlege
      real *8 wlege(*)
      real *8, allocatable :: w(:)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c

c     carve up workspace:
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
      lfr=(nterms+3)
c
      lused=ifr + lfr
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

      do i = 1, ns

      call l3dformmp0_charge_trunc(ier,rscale,sources(1,i),charge(i),
     1   center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),w(ifr))

      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp0_charge_trunc(ier,rscale,source,charge,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a charge located at the point SOURCE.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     charge  : charge strengths
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c     sss     : sign mutiplier array (see getsgnformpmp)
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
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 charge
      real *8 pp(0:nterms,0:nterms)
      real *8 fr(0:nterms+1)
      real *8 cscale,cscale1,cscale2,cscale3,rtmp
      real *8 a22,a21,a20,rmul,dd, a11,a10
      real *8 r,theta,phi,ctheta,stheta,cphi,sphi
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1),ephi1
      complex *16 eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      real *8 wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert dipole vector contributions to standard
c     n=0 moments about source position using standard
c     d+,d-,dz operators.
c
c     now shift n=1 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
      ephi(1) = dcmplx(cphi,sphi)
      ephi(-1) = dconjg(ephi(1))
c
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
c       
      call ylgndrufw(nterms,ctheta,pp,wlege,nlege)
c
        do i=0,nterms
        fr(i)=fr(i)*charge
        enddo
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
c       ... optimized for speed
c       use symmetries for real valued charges
c
            do ll = 0,nterms

               mpole(ll,0) = mpole(ll,0) + fr(ll)*pp(ll,0)

               do mm = 1,ll

               rtmp=fr(ll)*pp(ll,mm)
               mpole(ll, mm)= mpole(ll, mm) + rtmp*ephi(-mm)
               mpole(ll,-mm)= mpole(ll,-mm) + rtmp*ephi(+mm)

               enddo
            enddo
c
      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dformmp_dipole_trunc(ier,rscale,sources,dipvec,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS real-valued
c     dipole sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipvec(3,ns)    : dipole vector directions
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused, isss, lsss
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      real *8 center(3),sources(3,ns)
      real *8 dipvec(3,ns)
      real *8 rscale
      complex *16 mpole(0:nterms,-nterms:nterms)
      real *8 binom(0:120,0:2)
      integer nlege
      real *8 wlege(*)
      real *8, allocatable :: w(:)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c

c     carve up workspace:
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
      lfr=(nterms+3)
c
      isss=ifr + lfr
      lsss = 3*(2*nterms+1)
c
      lused=isss + lsss
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

ccc      call getsgnformpmp_dipole(w(isss),nterms)

      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*i)
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo


      do i = 1, ns

      call l3dformmp0_dipole_trunc(ier,rscale,sources(1,i),dipvec(1,i),
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,w(isss))

      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp0_dipole_trunc(ier,rscale,source,dipvec,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr,binom,sss)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a dipole located at the point SOURCE.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     dipvec  : dipole vector
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c     sss     : sign mutiplier array (see getsgnformpmp)
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
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      real *8 binom(0:120,0:2)
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 fr(0:nterms+1)
      real *8 cscale,cscale1,cscale2,cscale3,rtmp
      real *8 sss(-1:1,-nterms:nterms)
      real *8 a22,a21,a20,rmul,dd, a11,a10
      real *8 r,theta,phi,ctheta,stheta,cphi,sphi
      complex *16 mp1(-1:1)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1),ephi1
      complex *16 eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      real *8 wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert dipole vector contributions to standard
c     n=1 moments about source position using standard
c     d+,d-,dz operators.
c
      mp1(+1) = (-dipvec(1) + dipvec(2)*eye) /sqrt(2.0d0)
      mp1(0) = dipvec(3)
      mp1(-1) = dconjg(mp1(+1))
c
c     now shift n=1 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
      ephi(1) = dcmplx(cphi,sphi)
      ephi(-1) = dconjg(ephi(1))
c
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
c       
c      do l = 0,nterms
c      fr(l)=fr(l)/sqrt(2*l+1.0d0) *rscale
c      enddo
c      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
      do l = 0,nterms
      fr(l)=fr(l) *rscale
      enddo
      call ylgndrufw(nterms,ctheta,pp,wlege,nlege)
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
        if( 1 .eq. 2 ) then
c
c       ... reference code
c
         do m = -1,1
            do ll = 0,nterms-1
               cscale = fr(ll)
               lnew = 1+ll
c
               cscale2 = binom(lnew-m,1-m)*binom(lnew+m,1+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1            pp(ll,0)*cscale2*mp1(m)*cscale*sss(m,0)
               do mm = 1,ll
                  mnew1 = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew1,1-m)*binom(lnew+mnew1,1+m)
                  cscale2 = cscale2*cscale*sss(m,+mm)
                  cscale3 = binom(lnew-mnew2,1-m)*binom(lnew+mnew2,1+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
                  
                  mpole(lnew,mnew1) = mpole(lnew,mnew1)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp1(m)
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp1(m)

               enddo
            enddo
         enddo
         endif
c
c
        if( 2 .eq. 2 ) then
c
c       ... optimized for speed
c       use symmetries for real valued dipoles
c
            do ll = 0,nterms-1

               lnew = 1+ll
c
               cscale2 = binom(lnew,1)*binom(lnew,1)
               mpole(lnew,0) = mpole(lnew,0) +
     1            fr(ll)*pp(ll,0)*cscale2*mp1(0)

               rtmp = fr(ll)*pp(ll,0)
               cscale2 = binom(lnew-1,0)*binom(lnew+1,2)
               ztmp = rtmp*cscale2*mp1(1)
               mpole(lnew,+1) = mpole(lnew,+1)+ztmp
               mpole(lnew,-1) = mpole(lnew,-1)+dconjg(ztmp)

               do mm = 1,ll

               ctmp=fr(ll)*pp(ll,mm)*ephi(-mm)

               mnew1 = +mm
               mnew2 = -mm
               cscale2 = binom(lnew-mnew1,1)*binom(lnew+mnew1,1)
               ztmp = cscale2*ctmp*mp1(0)
               mpole(lnew,mnew1) = mpole(lnew,mnew1)+ztmp
               mpole(lnew,mnew2) = mpole(lnew,mnew2)+dconjg(ztmp)

               do m = 1,1
c
c       ... skip zero valued modes
c       
                  if( abs(mp1(m)) .eq. 0 ) cycle

                  mnew1 = m+mm
                  mnew2 = m-mm

ccc                  cscale2 = +binom(lnew-mnew1,0)*binom(lnew+mnew1,2)
ccc                  cscale3 = -binom(lnew-mnew2,0)*binom(lnew+mnew2,2)
                  cscale2 = +binom(lnew+mnew1,2)
                  cscale3 = -binom(lnew+mnew2,2)

                  ztmp0=cscale2*mp1(m)*ctmp
                  ztmp1=cscale3*mp1(m)*dconjg(ctmp)

                  mpole(lnew,+mnew1) = mpole(lnew,+mnew1)+ztmp0
                  mpole(lnew,+mnew2) = mpole(lnew,+mnew2)+ztmp1

                  mpole(lnew,-mnew2) = mpole(lnew,-mnew2)+dconjg(ztmp1)
                  mpole(lnew,-mnew1) = mpole(lnew,-mnew1)+dconjg(ztmp0)

               enddo
               enddo
            enddo
         endif
c
c
      return
      end
c
c
c
c
C***********************************************************************
      subroutine getsgnformpmp_dipole(sss,nterms)
C***********************************************************************
c
c     This subroutine creates a multiplier that holds the 
c     appropriate power of (-1) that arises in multipole-multipole shifts.
c     (See Greengard - 
c          Rapid Evaluation of Potential Fields... for notation).
c-------------------------------------------------------------------------
      implicit none
      integer m,mm,nterms
      real *8 sss(-1:1,-nterms:nterms)
c
      do m = -1,1
      do mm = -nterms,nterms
         sss(m,mm) = 1.0d0
      enddo
      enddo
      do m = -1,1
      do mm = -nterms,nterms
         if (( m .lt. 0) .and. (mm .gt. 0)) then
            if ( mm .le. -m ) sss(m,mm) = (-1)**mm
            if ( mm .gt. -m ) sss(m,mm) = (-1)**m
         endif
         if (( m .gt. 0) .and. (mm .lt. 0)) then
            if ( m .le. -mm ) sss(m,mm) = (-1)**m
            if ( m .gt. -mm ) sss(m,mm) = (-1)**mm
         endif
      enddo
      enddo
      return
      end



