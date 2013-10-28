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
c
c    l3dtaevalhessd  uses direct translation (not rotation/zshift)
c    and is reasonably optimized, precomputing the array of 
c    binomial/factorial terms that appear in the shift operator.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine l3dtaevalhessd_trunc(rscale,center,local,nterms,targ,
     1            pot,iffld,fld,ifhess,hess,scarray,wlege,nlege)
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the potential due to a local harmonic expansion.
c
c     pot =  sum sum  local(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     local  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     targ   :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c    scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                   L3DTAEVALHESSDINI(nterms,scarray))
c                   with dimension of scarray at least 10*(nterms+2)**2
c                   If nterms is changed, 
c                   l3dtaevalhessdini must be called again.
c   
c     OUTPUT:
c
c     pot    :    potential
c     fld    :    if (iffld .eq.1)   
c     hess   :    if (ifhess .eq.1)   
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c--------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ifhess,nlege
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,ldiff
      real *8 center(3),targ(3),wlege(nlege)
      real *8 zdiff(3)
      real *8 scarray(*),rscale
      real *8 cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 local2(0:2,-2:2)
      complex *16 z0,ima,ephi1,pot,fld(3),hess(6)
c
      real *8, allocatable :: pp(:,:)
      real *8, allocatable :: powers(:)
      complex *16, allocatable :: ppc(:,:)
      complex *16, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms,0:nterms))
      allocate(ppc(0:nterms,-nterms:nterms))
      allocate(powers(0:nterms+1))
      allocate(ephi(-nterms-1:nterms+1))
c
c     determine order of shifted expansion 
c
      ll = 0
      if (iffld.eq.1) ll = 1
      if (ifhess.eq.1) ll = 2
c
      do l = 0,ll
      do m = -l,l
         local2(l,m) = 0.0d0
      enddo
      enddo
c
      zdiff(1) = center(1) - targ(1)
      zdiff(2) = center(2) - targ(2)
      zdiff(3) = center(3) - targ(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
      dd = d*rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+1
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
      enddo
c
      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
      do l = 0,nterms
         do k = -l,l
            ppc(l,k) = pp(l,abs(k))*powers(l)*ephi(k)
         enddo
      enddo
c
c     shift to local expansion of order ll about target point
c
c     old code (not unrolled).
c
ccc      iuse = 1
ccc      do j = 0,ll
ccc         do k = -j,j
ccc            do l = j,nterms
ccc               do m = -l,l
ccc                  ldiff = l-j
ccc                  mm = m-k
ccc                  if (abs(mm).le.ldiff) then
ccc                     local2(j,k) = local2(j,k) +
ccc     1               ppc(ldiff,mm)*scarray(iuse)*local(l,m)
ccc                  endif
ccc                  iuse = iuse+1
ccc               enddo
ccc            enddo
ccc         enddo
ccc      enddo
c
      iuse = 1
      do l = 0,nterms
         do m = -l,l
            local2(0,0) = local2(0,0) +
     1      ppc(l,m)*scarray(iuse)*local(l,m)
            iuse = iuse+1
         enddo
      enddo
      if (ll.ge.1) then
         do l = 1,nterms
            ldiff = l-1
            do m = -l,l-2
               mm = m+1
               local2(1,-1) = local2(1,-1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               local2(1,0) = local2(1,0) +
     1         ppc(ldiff,m)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+2,l
               mm = m-1
               local2(1,1) = local2(1,1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
         enddo
      endif
      if (ll.eq.2) then
         do l = 2,nterms
            ldiff = l-2
            do m = -l,l-4
               mm = m+2
               local2(2,-2) = local2(2,-2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+1,l-3
               mm = m+1
               local2(2,-1) = local2(2,-1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               mm = m
               local2(2,0) = local2(2,0) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+3,l-1
               mm = m-1
               local2(2,1) = local2(2,1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+4,l
               mm = m-2
               local2(2,2) = local2(2,2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
         enddo
      endif
c
      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
      pot = local2(0,0)
c
c     fld comes from l=1 modes
c
      rfac = rscale/sqrt(2.0d0)
      if (iffld.eq.1) then
         fld(3) = -local2(1,0)*rscale
         fld(1) = rfac*(local2(1,1) + local2(1,-1))
         fld(2) = rfac*ima*(local2(1,1) - local2(1,-1))
      endif
c
c     hess comes from l=2 modes
c
      if (ifhess.eq.1) then
ccc         rfac = rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
         rfac = rfac*rscale*sqrt(3.0d0)
         z0 = rscale*rscale*local2(2,0)
         hess(1) = rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(2) = -rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(3) = 2*z0
         hess(4) = rfac*ima*(local2(2,2) - local2(2,-2))
         hess(5) = -rfac*(local2(2,1) + local2(2,-1))
         hess(6) = -rfac*ima*(local2(2,1) - local2(2,-1))
      endif
c
      return
      end
c
c
c
c
c
