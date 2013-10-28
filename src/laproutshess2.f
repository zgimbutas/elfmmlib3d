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
c
c
c
c
c**********************************************************************
      subroutine l3dmpevalhess(rscale,center,mpole,nterms,ztarg,
     1		pot,iffld,fld,ifhess,hess,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
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
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     hess   :    Hessian at ztarg (if requested)
c                 in order (dxx,dyy,dzz,dxy,dxz,dyz)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms,iffld,ifhess
      integer nterms2,iloc,lloc,lused
      real *8 rscale,center(3),ztarg(3)
      real *8, allocatable :: w(:)
      complex *16 pot,fld(3),hess(6)
      complex *16 mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
c     for local expansion
c
      nterms2=0
      if (iffld.eq.1) nterms2=1
      if (ifhess.eq.1) nterms2=2
      iloc=1
      lloc=2*(nterms2+1)*(2*nterms2+1)+5
c
      lused=iloc+lloc
      allocate(w(lused))
c
      call l3dmplocquadu(rscale,center,mpole,nterms,rscale,
     1	   ztarg,w(iloc),nterms2,ier)
c
      call l3devalhessloc(rscale,w(iloc),nterms2,
     $   pot,iffld,fld,ifhess,hess)
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
      subroutine l3dtaevalhess(rscale,center,mpole,nterms,ztarg,
     1		pot,iffld,fld,ifhess,hess,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the 
c     potential due to a local multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  r^{n}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
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
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     hess   :    Hessian at ztarg (if requested)
c                 in order (dxx,dyy,dzz,dxy,dxz,dyz)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms,iffld,ifhess
      integer nterms2,iloc,lloc,lused
      real *8 rscale,center(3),ztarg(3)
      real *8, allocatable :: w(:)
      complex *16 pot,fld(3),hess(6)
      complex *16 mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
c     for local expansion
c
      nterms2=0
      if (iffld.eq.1) nterms2=1
      if (ifhess.eq.1) nterms2=2
      iloc=1
      lloc=2*(nterms2+1)*(2*nterms2+1)+5
c
      lused=iloc+lloc
      allocate(w(lused))
c
      call l3dloclocquadu(rscale,center,mpole,nterms,rscale,
     1	   ztarg,w(iloc),nterms2,ier)
c
      call l3devalhessloc(rscale,w(iloc),nterms2,
     $   pot,iffld,fld,ifhess,hess)
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine l3devalhessloc(rscale,local,nterms,
     $     pot,iffld,fld,ifhess,hess)
c**********************************************************************
c
c     Locally u = local(0,0) +
c                 local(1,-1) * r  * Y_{1,-1} +
c                 local(1,0)  * r  * Y_{1,0}  +
c                 local(1,1)  * r  * Y_{1,1}  +
c                 local(2,-2) * r2 * Y_{2,-2} +
c                 local(2,-1) * r2 * Y_{2,-1} +
c                 local(2,0)  * r2 * Y_{2,0}  +
c                 local(2,1)  * r2 * Y_{2,1}  +
c                 local(2,2)  * r2 * Y_{2,2}  +
c     so evaluation is easy (depending only on scaling convention
c     for definitions of Y_{n,m}. 
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ifhess
      real *8 rscale,pi,rfac
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 pot,fld(3),hess(6),eye,z0
c
      eye = dcmplx(0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
      pot = local(0,0)
c
c     fld comes from n=1 modes
c
      if( iffld .eq. 1 ) then
      rfac = sqrt(2.0d0)*rscale
      fld(1) = rfac*(local(1,1) + local(1,-1))/2.0d0
      fld(2) = rfac*eye*(local(1,1) - local(1,-1))/2.0d0
      fld(3) = -local(1,0)*rscale
      endif
c
c     hess comes from n=2 modes
c
      if( ifhess .eq. 1 ) then
      rfac = rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
      z0 = rscale*rscale*local(2,0)
      hess(1) = rfac*(local(2,2) + local(2,-2)) - z0
      hess(2) = -rfac*(local(2,2) + local(2,-2)) - z0
      hess(3) = 2*z0
      hess(4) = rfac*eye*(local(2,2) - local(2,-2))
      hess(5) = -rfac*(local(2,1) + local(2,-1))
      hess(6) = -rfac*eye*(local(2,1) - local(2,-1))
      endif
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dallhess(iffld,ifhess,sources,charge,ns,
     1                   target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field FLD,
c     and Hessian HESS at the target point TARGET, due to a collection 
c     of charges at SOURCE(3,ns). 
c     
c              	pot =  sum 1/r
c		fld =  -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> dont compute 
c		                   ifhess = 1 -> do compute 
c     sources(3,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)        : calculated potential
c     fld   (complex *16)        : calculated gradient
c     hess  (complex *16)        : calculated Hessian
c
c---------------------------------------------------------------------
      implicit none
      integer iffld,ifhess,ns,i,j
      real *8 sources(3,ns),target(3)
      complex *16 pot,fld(3),hess(6),potloc,fldloc(3),hessloc(6)
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
      if (ifhess.eq.1) then
         do i = 1,6
            hess(i) = 0.0d0
         enddo
      endif
c
      do i = 1,ns
         call lpotfld3dhess(iffld,ifhess,sources(1,i),charge(i),
     1        target,potloc,fldloc,hessloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
         if (ifhess.eq.1) then
            do j = 1,6
               hess(j) = hess(j) + hessloc(j)
            enddo
         endif
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotfld3dhess(iffld,ifhess,source,charge,target,
     1                         pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field FLD
c     and Hesian HESS at the target point TARGET, due to a charge at 
c     SOURCE. 
c     
c              	pot  = 1/r
c		fld  = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld     : flag for computing gradient
c	                 	iffld = 0 -> dont compute 
c		                iffld = 1 -> do compute 
c     ifhess    : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     fld       : calculated gradient
c     hess      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit none
      integer iffld,ifhess
      real *8 source(3),target(3)
      real *8 xdiff,ydiff,zdiff,dd,d,dinv,ddinv,dddinv,dddddinv
      complex *16 pot,fld(3),hess(6)
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
c
c       O(1/r)
c
      pot=charge*dinv
c
      ddinv=dinv*dinv
      dddinv=dinv*ddinv
      dddddinv=ddinv*dddinv
c
      if (iffld.eq.1) then
c
c       O(1/r^2)
c
         cd=charge*dddinv
         fld(1)=cd*xdiff
         fld(2)=cd*ydiff
         fld(3)=cd*zdiff
      endif
c
      if (ifhess.eq.1) then
c
c       O(1/r^3)
c
         cd=charge*dddddinv
         hess(1)=cd*(3*xdiff*xdiff-dd)
         hess(2)=cd*(3*ydiff*ydiff-dd)
         hess(3)=cd*(3*zdiff*zdiff-dd)
         hess(4)=cd*3*xdiff*ydiff
         hess(5)=cd*3*xdiff*zdiff
         hess(6)=cd*3*ydiff*zdiff
      endif
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dallhess_dp(iffld,ifhess,
     $     sources,dipstr,dipvec,ns,target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT field FLD and
c     Hessian HESS at the target point TARGET, due to a collection 
c     of dipoles at SOURCE(3,ns). 
c     
c              	pot = (dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c		fld = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     ifhess       : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
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
c     hess         : calculated hessian
c
c----------------------------------------------------------------------
      implicit none
      integer iffld,ifhess,ns,i,j
      real *8 sources(3,ns),target(3)
      real *8 dipvec(3,ns)
      complex *16 pot,fld(3),hess(6),potloc,fldloc(3),hessloc(6)
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
      if (ifhess.eq.1) then
         do i = 1,6
            hess(i) = 0.0d0
         enddo
      endif
c
      do i = 1,ns
         call lpotfld3dhess_dp(iffld,ifhess,
     $        sources(1,i),dipstr(i),dipvec(1,i),
     1        target,potloc,fldloc,hessloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
         if (ifhess.eq.1) then
            do j = 1,6
               hess(j) = hess(j) + hessloc(j)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dhess_dp(
     $     iffld,ifhess,source,dipstr,dipvec,target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT field FLD and
c     Hessian HESS at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = (dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c		fld = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     ifhess       : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
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
c     hess         : calculated hessian
c
c----------------------------------------------------------------------
      implicit none
      integer iffld,ifhess
      real *8 source(3),target(3)
      real *8 dipvec(3)
      real *8 xdiff,ydiff,zdiff,dd,d,dinv,ddinv,dddinv,dotprod
      real *8 rx,ry,rz,dx,dy,dz,rtmp,dddddinv
      complex *16 pot,fld(3),hess(6)
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
      ddinv=dinv*dinv  
      dddinv=ddinv*dinv  
c
c       O(1/r^2)
c
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)+zdiff*dipvec(3)
      cd = dipstr*dotprod
      pot=cd*dddinv
c
      dddddinv=ddinv*dddinv
c
      if (iffld.eq.1) then
c
c       O(1/r^3)
c
         ztttt = 3.0d0*cd*dddddinv
         cd2 = dipstr*dddinv
         fld(1)=ztttt*xdiff-cd2*dipvec(1)
         fld(2)=ztttt*ydiff-cd2*dipvec(2)
         fld(3)=ztttt*zdiff-cd2*dipvec(3)
      endif 
c
      if (ifhess.eq.1) then
c
c       O(1/r^4)
c
         rx=xdiff
         ry=ydiff
         rz=zdiff
c
         dx=rx*dinv
         dy=ry*dinv
         dz=rz*dinv
c
         rtmp=dotprod
c
         hess(1)=3*(rtmp*(5*dx*dx-1)-(dipvec(1)*rx+dipvec(1)*rx))
         hess(2)=3*(rtmp*(5*dy*dy-1)-(dipvec(2)*ry+dipvec(2)*ry))
         hess(3)=3*(rtmp*(5*dz*dz-1)-(dipvec(3)*rz+dipvec(3)*rz))
c
         hess(4)=3*(rtmp*(5*dx*dy)-(dipvec(2)*rx+dipvec(1)*ry))
         hess(5)=3*(rtmp*(5*dx*dz)-(dipvec(3)*rx+dipvec(1)*rz))
         hess(6)=3*(rtmp*(5*dy*dz)-(dipvec(3)*ry+dipvec(2)*rz))
c
         cd=dipstr*dddddinv
         hess(1)=hess(1)*cd
         hess(2)=hess(2)*cd
         hess(3)=hess(3)*cd
         hess(4)=hess(4)*cd
         hess(5)=hess(5)*cd
         hess(6)=hess(6)*cd
c
      endif 
c
      return
      end
cc
c
c
c
