cc Copyright (C) 2011: Leslie Greengard and Zydrunas Gimbutas
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
c      forming multipole expansions due to quadrupole
c      sources.
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
c
c-----------------------------------------------------------------------
c
c      f90 version, using allocate
c
c      L3DFORMMP_QUAD: creates multipole expansion (outgoing) due to 
c                 a collection of quadrupoles.
c                 (calls L3DFORMMP1/L3DFORMMP0_QUAD )
c
c      LPOTFLD3DALL_QUAD: 
c                 direct calculation for a collection of quadrupoles
c      LPOTFLD3D_QUAD : direct calculation for a single quadrupole
c
c
c      L3DFORMTA_QUAD: NOT YET IMPLEMENTED
c                 (calls L3DFORMTA1/L3DFORMTA0_QUAD )
c
c
C***********************************************************************
      subroutine l3dformmp_quad(ier,rscale,sources,quadvec,ns,
     1                  center,nterms,mpole)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     quadrupoles sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     quadvec(6,ns)    : quadrupoles vector direction 
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
      integer nterms,ns,i,l,m, ier, lused, istotal
      real *8 center(3),sources(3,ns)
      real *8 quadvec(6,ns)
      real *8 rscale
      real *8, allocatable :: sss(:)
      complex *16 mpole(0:nterms,-nterms:nterms)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      istotal = 5*(2*nterms+1)
      allocate(sss(istotal))
c
      call getsgnformpmp_quad(sss,nterms)
      do i = 1, ns
         call l3dformmp1_quad(ier,rscale,sources(1,i),
     1        quadvec(1,i),center,nterms,mpole,sss)
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_quad(ier,rscale,source,quadvec,
     1		center,nterms,mpole,sss)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a quadrupole located at the point SOURCE.
c     This is a memory management routine. Work is done in the
c     secondary call to l3dformmp0_quad below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     quadvec  : quadrupole vector
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c     sss     : sign mutiplier array (see getsgnformpmp_quad)
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
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr,lused
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 quadvec(6)
      real *8 sss(-2:2,-nterms:nterms)
      complex *16 mpole(0:nterms,-nterms:nterms)
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
c
      call l3dformmp0_quad(jer,rscale,source,quadvec,
     1		center,nterms,mpole,w(ipp),w(iephi),
     2          w(ifr),sss)
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_quad(ier,rscale,source,quadvec,
     1		center,nterms,mpole,pp,ephi,fr,sss)
c**********************************************************************
c
c     See l3dformmp1_quad for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew,mnew2
      real *8 binom(0:120,0:4)
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 quadvec(6)
      real *8 pp(0:nterms,0:nterms)
      real *8 fr(0:nterms+1)
      real *8 cscale,cscale2,cscale3
      real *8 sss(-2:2,-nterms:nterms)
      real *8 a22,a21,a20,rmul,dd
      real *8 r,theta,phi,ctheta,stheta,cphi,sphi
      complex *16 mp2(-2:2)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1),ephi1
      complex *16 eye,zmul
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
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
      ephi(-1) = dcmplx(cphi,-sphi)
c
      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*(i))
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo
      do i = 3,2*nterms
         binom(i,3) = dsqrt((i*(i-1)*(i-2))/6.0d0)
      enddo
      do i = 4,2*nterms
         binom(i,4) = dsqrt((i*(i-1)*(i-2)*(i-3))/24.0d0)
      enddo
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
C
      call ylgndr(nterms,ctheta,pp)
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
         do m = -2,2
            do ll = 0,nterms-2
               cscale = rscale*rscale*fr(LL)/dsqrt(2*ll+1.0d0)
               lnew = 2+ll
c
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1         pp(ll,0)*cscale2*mp2(m)*cscale*sss(m,0)
               do mm = 1,ll
                  mnew = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew,2-m)*binom(lnew+mnew,2+m)
                  cscale2 = cscale2*cscale*sss(m,mm)
                  cscale3 = binom(lnew-mnew2,2-m)*binom(lnew+mnew2,2+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
c
                  mpole(lnew,mnew) = mpole(lnew,mnew)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp2(m)
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp2(m)
               enddo
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
      subroutine l3dformmp_quad_trunc(ier,rscale,sources,quadvec,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     quadrupoles sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     quadvec(6,ns)    : quadrupoles vector direction 
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
      real *8 quadvec(6,ns)
      real *8 rscale
      complex *16 mpole(0:nterms,-nterms:nterms)
      real *8 binom(0:120,0:4)
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
      lsss = 5*(2*nterms+1)
c
      lused=isss + lsss
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

ccc      call getsgnformpmp_quad(w(isss),nterms)

      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*i)
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo
      do i = 3,2*nterms
         binom(i,3) = dsqrt((i*(i-1)*(i-2))/6.0d0)
      enddo
      do i = 4,2*nterms
         binom(i,4) = dsqrt((i*(i-1)*(i-2)*(i-3))/24.0d0)
      enddo


      do i = 1, ns

c        call l3dformmp1_quad(ier,rscale,sources(1,i),
c     1        quadvec(1,i),center,nterms,mpole,w(isss))

      call l3dformmp0_quad_trunc(ier,rscale,sources(1,i),quadvec(1,i),
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,w(isss))

      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_quad_trunc(ier,rscale,source,quadvec,
     1		center,nterms,mpole,wlege,nlege,binom,sss)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a quadrupole located at the point SOURCE.
c     This is a memory management routine. Work is done in the
c     secondary call to l3dformmp0_quad below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     quadvec  : quadrupole vector
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
      integer ier,nterms
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr,lused
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 quadvec(6)
      real *8 binom(0:120,0:4)
      real *8 sss(-2:2,-nterms:nterms)
      complex *16 mpole(0:nterms,-nterms:nterms)
      integer nlege
      real *8 wlege(*)
c
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

      call l3dformmp0_quad_trunc(jer,rscale,source,quadvec,
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,sss)

      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_quad_trunc(ier,rscale,source,quadvec,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr,binom,sss)
c**********************************************************************
c
c     See l3dformmp1_quad for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      real *8 binom(0:120,0:4)
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 quadvec(6)
      real *8 pp(0:nterms,0:nterms)
      real *8 fr(0:nterms+1)
      real *8 cscale,cscale1,cscale2,cscale3,rtmp
      real *8 sss(-2:2,-nterms:nterms)
      real *8 a22,a21,a20,rmul,dd
      real *8 r,theta,phi,ctheta,stheta,cphi,sphi
      complex *16 mp2(-2:2)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1),ephi1
      complex *16 eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      real *8 wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
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
c      fr(l)=fr(l)/sqrt(2*l+1.0d0) *rscale*rscale
c      enddo
c      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
      do l = 0,nterms
      fr(l)=fr(l) *rscale*rscale
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
         do m = -2,2
            do ll = 0,nterms-2
               cscale = fr(ll)
               lnew = 2+ll
c
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1            pp(ll,0)*cscale2*mp2(m)*cscale*sss(m,0)
               do mm = 1,ll
                  mnew1 = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew1,2-m)*binom(lnew+mnew1,2+m)
                  cscale2 = cscale2*cscale*sss(m,+mm)
                  cscale3 = binom(lnew-mnew2,2-m)*binom(lnew+mnew2,2+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
                  
                  mpole(lnew,mnew1) = mpole(lnew,mnew1)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp2(m)
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp2(m)

               enddo
            enddo
         enddo
         endif
c
c
        if( 2 .eq. 2 ) then
c
c       ... optimized for speed
c       use symmetries for real valued quadrupoles
c
            do ll = 0,nterms-2

               lnew = 2+ll
c
               cscale2 = binom(lnew,2)*binom(lnew,2)
               mpole(lnew,0) = mpole(lnew,0) +
     1            fr(ll)*pp(ll,0)*cscale2*mp2(0)

               rtmp = fr(ll)*pp(ll,0)
               do m = 1,2
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               ztmp = rtmp*cscale2*mp2(m)
               mpole(lnew,+m) = mpole(lnew,+m)+ztmp
               mpole(lnew,-m) = mpole(lnew,-m)+dconjg(ztmp)
               enddo

               do mm = 1,ll

               ctmp=fr(ll)*pp(ll,mm)*ephi(-mm)

               mnew1 = +mm
               mnew2 = -mm
               cscale2 = binom(lnew-mnew1,2)*binom(lnew+mnew1,2)
               ztmp = cscale2*ctmp*mp2(0)
               mpole(lnew,mnew1) = mpole(lnew,mnew1)+ztmp
               mpole(lnew,mnew2) = mpole(lnew,mnew2)+dconjg(ztmp)

               do m = 1,2
c
c       ... skip zero valued modes
c       
                  if( abs(mp2(m)) .eq. 0 ) cycle

                  mnew1 = m+mm
                  mnew2 = m-mm

                  if( m .eq. 1 ) then
                  cscale2 = +binom(lnew-mnew1,1)*binom(lnew+mnew1,3)
                  cscale3 = -binom(lnew-mnew2,1)*binom(lnew+mnew2,3)
                  else
ccc                  cscale2 = +binom(lnew-mnew1,0)*binom(lnew+mnew1,4)
ccc                  cscale3 = +binom(lnew-mnew2,0)*binom(lnew+mnew2,4)
                  if( mm .eq. 1 ) then
                  cscale2 = +binom(lnew+mnew1,4)
                  cscale3 = -binom(lnew+mnew2,4)
                  else
                  cscale2 = +binom(lnew+mnew1,4)
                  cscale3 = +binom(lnew+mnew2,4)
                  endif
                  endif

                  ztmp0=cscale2*mp2(m)*ctmp
                  ztmp1=cscale3*mp2(m)*dconjg(ctmp)

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
c**********************************************************************
      subroutine lpotfld3dall_quad(iffld,sources,quadvec,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of quadrupoles 
c     at SOURCE(3,ns). 
c     
c		fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     quadvec(6,ns)  : dipole direction
c     ns            : number of sources
c     charge(ns)    : charge strength
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,i,ns
      real *8 sources(3,ns),target(3)
      real *8 quadvec(6,ns)
      complex *16 pot,fld(3),potloc,fldloc(3)
c
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotfld3d_quad(iffld,sources(1,i),quadvec(1,i),
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
      subroutine lpotfld3d_quad(iffld,source,quadvec,target,
     1                        pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a quadrupole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c    
c               pot = quadvec(1)*V_xx +
c                     quadvec(2)*V_yy +
c                     quadvec(3)*V_zz +
c                     quadvec(4)*V_xy +
c                     quadvec(5)*V_xz +
c                     quadvec(6)*V_yz 
c
c      V_xx = (-1/r^3 + 3*dx**2/r^5)
c      V_xy = 3*dx*dy/r^5
c      V_xz = 3*dx*dz/r^5
c      V_yy = (-1/r^3 + 3*dy**2/r^5)
c      V_yz = 3*dy*dz/r^5
c      V_zz = (-1/r^3 + 3*dz**2/r^5)
c
c		fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     source(3)    : location of the source 
c     quadvec(3,ns) : quadrupole vector
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
      real *8 quadvec(6),rr(3)
      real *8 cd,d5,d3,v11,v12,v13,v22,v23,v33
      real *8 d7,v111,v112,v113,v122,v123,v133
      real *8 v222,v223,v233,v333
      complex *16 pot,fld(3)
c
c
c ... Caculate offsets and distance
c
      rr(1)=target(1)-source(1)
      rr(2)=target(2)-source(2)
      rr(3)=target(3)-source(3)
      cd = sqrt(rr(1)*rr(1)+rr(2)*rr(2)+rr(3)*rr(3))
      d5 = 1.0d0/(cd**5)
      d3 = 1.0d0/(cd**3)
c
c ... Calculate the potential and field in the regular case:
c
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      v11 = -d3 + 3*d5*rr(1)**2
      v12 = 3*d5*rr(1)*rr(2)
      v13 = 3*d5*rr(1)*rr(3)
      v22 = -d3 + 3*d5*rr(2)**2
      v23 = 3*d5*rr(2)*rr(3)
      v33 = -d3 + 3*d5*rr(3)**2
c
      pot=v11*quadvec(1)+v22*quadvec(2)+v33*quadvec(3)
      pot=pot+v12*quadvec(4)+v13*quadvec(5)+v23*quadvec(6)
ccc      call prin2(' cd is *',cd,1)
ccc      call prin2(' d3 is *',d3,1)
ccc      call prin2(' v11 is *',v11,1)
ccc      call prin2(' v22 is *',v22,1)
ccc      call prin2(' v33 is *',v33,1)
ccc      call prin2(' v12 is *',v12,1)
ccc      call prin2(' v13 is *',v13,1)
ccc      call prin2(' v23 is *',v23,1)
ccc      call prin2(' pot is *',pot,1)
      if (iffld.eq.1) then
         d7 = 1.0d0/(cd**7)
         v111 = 9*rr(1)*d5 - 15*d7*rr(1)**3
         v112 = 3*rr(2)*d5 - 15*d7*rr(2)*rr(1)**2
         v113 = 3*rr(3)*d5 - 15*d7*rr(3)*rr(1)**2
         v122 = 3*rr(1)*d5 - 15*d7*rr(1)*rr(2)**2
         v123 = -15*d7*rr(1)*rr(2)*rr(3)
         v133 = 3*rr(1)*d5 - 15*d7*rr(1)*rr(3)**2
         v222 = 9*rr(2)*d5 - 15*d7*rr(2)**3
         v223 = 3*rr(3)*d5 - 15*d7*rr(3)*rr(2)**2
         v233 = 3*rr(2)*d5 - 15*d7*rr(2)*rr(3)**2
         v333 = 9*rr(3)*d5 - 15*d7*rr(3)**3
ccc         call prin2(' v111 is *',v111,1)
ccc         call prin2(' v112 is *',v112,1)
ccc         call prin2(' v113 is *',v113,1)
ccc         call prin2(' v122 is *',v122,1)
ccc         call prin2(' v123 is *',v123,1)
ccc         call prin2(' v133 is *',v133,1)
ccc         call prin2(' v222 is *',v222,1)
ccc         call prin2(' v223 is *',v223,1)
ccc         call prin2(' v233 is *',v233,1)
ccc         call prin2(' v333 is *',v333,1)
         fld(1)=v111*quadvec(1)+v122*quadvec(2)+v133*quadvec(3)
         fld(1)=fld(1)+v112*quadvec(4)+v113*quadvec(5)+v123*quadvec(6)
         fld(2)=v112*quadvec(1)+v222*quadvec(2)+v233*quadvec(3)
         fld(2)=fld(2)+v122*quadvec(4)+v123*quadvec(5)+v223*quadvec(6)
         fld(3)=v113*quadvec(1)+v223*quadvec(2)+v333*quadvec(3)
         fld(3)=fld(3)+v123*quadvec(4)+v133*quadvec(5)+v233*quadvec(6)
         fld(1) = -fld(1)
         fld(2) = -fld(2)
         fld(3) = -fld(3)
      endif 
      return
      end
cc
c
c
c**********************************************************************
      subroutine l3dformta_quad(ier,rscale,sources,quadvec,ns,
     1		           center,nterms,locexp)
c**********************************************************************
c
c    NOT YET IMPLEMENTED
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS quadrupoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1_quad/l3dformta0_quad below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     quadvec    : quadrupole vector
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c     w         : workspace
c     lw        : workspace length
c                           at least (nterms+1)**2 +
c                              6*(nterms+1) + 7000
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c
c     locexp    : coeffs for the j-expansion
c     lused     : amount of work space "w" used
c
c
c----------------------------------------------------------------------
      implicit none
      integer ier,ns,nterms,i,l,m
      real *8 rscale,sources(3,ns),center(3),rs
      real *8 quadvec(6,ns)
      complex *16 locexp(0:nterms,-nterms:nterms)
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
         call l3dformta1_quad(ier,rscale,sources(1,i),
     1		quadvec(1,i),center,nterms,locexp)
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
      subroutine l3dformta1_quad(ier,rscale,source,quadvec,
     &		center,nterms,locexp)
c**********************************************************************
c
c    NOT YET IMPLEMENTED
c     This subroutine creates the local expansion about CENTER
c     due to a single quadrupole located at SOURCE.
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
c     quadvec   : quadrupole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c     w         : workspace
c     lw        : workspace length
c                           at least (nterms+1)**2 +
c                                    6*(nterms+1) + 1000
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   ier=2 insufficient memory in workspace w
c     locexp : coefficients of the local expansion
c     lused  : amount of work space "w" used
c
c
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms,ipp,lpp,iephi,lephi,ifr,lfr,lused
      real *8 rscale,source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 quadvec(6)
      complex *16 locexp(0:nterms,-nterms:nterms)
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
      call l3dformta0_quad(rscale,source,quadvec,
     &   center,nterms,locexp,w(ipp),w(iephi),w(ifr))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0_quad(rscale,source,quadvec,
     &		center,nterms,locexp,pp,ephi,fr)
c**********************************************************************
c
c    NOT YET IMPLEMENTED
c     See l3dformta_quad/l3dformta1_quad for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,i
      real *8 rscale,source(3),center(3),zdiff(3)
      real *8 quadvec(6)
      real *8 pp(0:nterms,0:nterms)
      real *8 done,r,theta,phi,ctheta,stheta,cphi,sphi,d
      complex *16 locexp(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fr(0:nterms+1),ztmp,z
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
ccc      do i=0,nterms
ccc         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
ccc      enddo
c
c    NOT YET IMPLEMENTED
c
      return
      end
c
c
C***********************************************************************
      subroutine getsgnformpmp_quad(sss,nterms)
C***********************************************************************
c
c     This subroutine creates a multiplier that holds the 
c     appropriate power of (-1) that arises in multipole-multipole shifts.
c     (See Greengard - 
c          Rapid Evaluation of Potential Fields... for notation).
c-------------------------------------------------------------------------
      implicit none
      integer m,mm,nterms
      real *8 sss(-2:2,-nterms:nterms)
c
      do m = -2,2
      do mm = -nterms,nterms
         sss(m,mm) = 1.0d0
      enddo
      enddo
      do m = -2,2
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
