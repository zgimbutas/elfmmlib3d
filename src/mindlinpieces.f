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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of Mindlin B and C parts for half space
c       elastostatic Green's functions in R^3.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c***********************************************************************
      subroutine intker_mindlinb(icomp,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,ztrg,pot,iffld,fld)
c***********************************************************************
c
c     This subroutine computes the displacement due to the 
c     harmonic part of the Mindlin B solution and all of its
c     derivatives, from which stress and strain can be extracted.
c     (See okada_imageanalysis for details.)
c     It doesn't include the 1/rmu scaling.
c
c     INPUT:
c
c     icomp      desired component of displacement vector
c     rlame(2)   Lame coefficients supplied in the form
c                rlame(1) = rlam, rlame(2) = rmu
c     source(3)  image source location.
c                This source must be the reflected image in the 
c                upper half-space, NOT the original source in the 
c                lower half-space.
c     ifsingle   (the single layer kernel flag)
c                0 means ignore the <<force>> vector 
c                1 means a <<force>> vector is supplied 
c     sigma_sl   force vector
c     ifdouble   (the double layer kernel flag)
c                0 means ignore the <<dislocation>> 
c                1 means a <<dislocation vector>> is supplied 
c                  with surface normal <<dipvec>>
c     dipstr     dislocation vector
c     dipvec     surface orientation
c                
c     ztrg       target location
c     pot        component <<icomp>> of displacement vector
c     iffld      fld computation flag: 
c                0 means don't compute
c                1 means compute fld = -gradient
c
c     fld        -gradient(pot) at ztrg
c-----------------------------------------------------------------------
      implicit none
      integer icomp,ifsingle,ifdouble,iffld,nder(3)
      integer j,k,l,n
      real *8 rlame(2),source(3),ztrg(3),sigma_sl(3)
      real *8 uder(3,3,3),uuu,uuuder(3)
      real *8 udertot(0:4,0:4,0:4)
      real *8 dipstr(3),dipvec(3),rr(3)
      real *8 denom1,denom2,rr3,cd,rlam,rmu,scale,sgn
      real *8 u11,u12,u13,u22,u23,u33,uder1,uder2,uder4
      real *8 alpha,d1,d2,dotprod,dx,dy,dz
      complex *16 pot,fld(3)
c      
c    
c
      rlam = rlame(1)
      rmu = rlame(2)
      alpha = (rlam+rmu)/(rlam+2*rmu)
      scale = (1.0d0-alpha)/alpha
c
      pot = 0.0d0
      fld(1) = 0.0d0
      fld(2) = 0.0d0
      fld(3) = 0.0d0
c
      dx = ztrg(1)-source(1)
      dy = ztrg(2)-source(2)
      dz = ztrg(3)-source(3)
      rr(1) = dx
      rr(2) = dy
      rr(3) = -dz
      rr3 = rr(3)
      cd = sqrt(dx*dx+dy*dy+dz*dz)
c
      if ( (ifdouble.eq.1) .or. (iffld.eq.1) ) then
c
c     Compute three derivatives of harmonic function 
c     U = R_3 log(R+R_3) - R
c     with respect to R_1, R_2 and R_3. We will use these
c     derivatives to compute both the gradient of the SLP
c     and the DLP.
c
c
         denom1 = ((cd-dz)*cd)**3
         denom2 = ((cd-dz)**2)*cd**3
         uder(1,1,1) = 3*dx*cd*((cd-dz)*cd-dx**2) + dz*dx**3
         uder(1,1,1) = uder(1,1,1)/denom1
         uder(1,1,2) = (cd-dz)*dy*(cd*cd-dx*dx)-2*cd*dy*dx**2
         uder(1,1,2) = uder(1,1,2)/denom1
         uder(1,1,3) = (cd-dz)*(cd*cd-dx*dx) - dx*dx*cd
         uder(1,1,3) = uder(1,1,3)/denom2
         uder(1,2,1) = uder(1,1,2)       
         uder(1,2,2) = (cd-dz)*dx*(cd*cd-dy*dy)-2*cd*dx*dy**2
         uder(1,2,2) = uder(1,2,2)/denom1
         uder(1,2,3) = -dx*dy*(2*cd-dz)
         uder(1,2,3) = uder(1,2,3)/denom2
         uder(1,3,1) = uder(1,1,3)       
         uder(1,3,2) = uder(1,2,3)       
         uder(1,3,3) = -dx*(cd-dz)**2
         uder(1,3,3) = uder(1,3,3)/denom2
c
         if (icomp.ne.1) then
            uder(2,1,1) = uder(1,1,2)
            uder(2,1,2) = uder(1,2,2)
            uder(2,1,3) = uder(1,2,3)
            uder(2,2,1) = uder(1,2,2)
            uder(2,2,2) = 3*dy*cd*((cd-dz)*cd-dy**2) + dz*dy**3
            uder(2,2,2) = uder(2,2,2)/denom1
            uder(2,2,3) = (cd-dz)*(cd*cd-dy*dy) - dy*dy*cd
            uder(2,2,3) = uder(2,2,3)/denom2
            uder(2,3,1) = uder(1,2,3)
            uder(2,3,2) = uder(2,2,3)
            uder(2,3,3) = -dy*(cd-dz)**2
            uder(2,3,3) = uder(2,3,3)/denom2
         endif
c
         if (icomp.eq.3) then
            uder(3,1,1) = uder(1,1,3)
            uder(3,1,2) = uder(1,2,3)
            uder(3,1,3) = uder(1,3,3)
            uder(3,2,1) = uder(1,2,3)
            uder(3,2,2) = uder(2,2,3)
            uder(3,2,3) = uder(2,3,3)
            uder(3,3,1) = uder(1,3,3)
            uder(3,3,2) = uder(2,3,3)
            uder(3,3,3) = dz/(cd**3)
         endif
      endif
c
      if (ifsingle.eq.1) then
c
c
c     Compute two derivatives of harmonic function 
c     U = R_3 log(R+R_3) - R
c     with respect to R_1, R_2 and R_3. These play a role only 
c     in the SLP.
c
         d1 = 1.0d0/(cd+rr(3))
         d2 = 1.0d0/(cd*(cd+rr(3)))
         u11 = d1*(-1.0d0 + d2*rr(1)**2)
         u12 = d1*d2*rr(1)*rr(2)
         u13 = d2*rr(1)
         u22 = d1*(-1.0d0 + d2*rr(2)**2)
         u23 = d2*rr(2)
         u33 = 1.0d0/cd
c
c     See notes (okada_imageanalysis) for explanation.
c
         if (icomp.eq.1) then
            pot = -sigma_sl(1)*u11 - sigma_sl(2)*u12 + sigma_sl(3)*u13
            pot = pot*scale
         else if (icomp.eq.2) then
            pot = -sigma_sl(1)*u12 - sigma_sl(2)*u22 + sigma_sl(3)*u23
            pot = pot*scale
         else 
            pot = -sigma_sl(1)*u13 - sigma_sl(2)*u23 + sigma_sl(3)*u33
            pot = pot*scale
         endif 
c
c      compute contribution to gradient of SLP for component icomp
c
c      Because Laplace FMM is based on FIELD rather than GRADIENT,
c      we compute FIELD = -GRADIENT here. We just take (-gradient)
c      of preceding formula. Since
c      d/dx_1 = d/dR_1, d/dx_2 = d/dR_2 and d/dx_3 = -d/dR_3
c      fld(1) and fld(2) and fld(3) use different signs below.
c
         if (iffld.eq.1) then
            sgn = scale
            do k = 1,3
               if (k.eq.3) sgn = -sgn
               fld(k) = fld(k) + sgn*(
     1          uder(icomp,1,k)*sigma_sl(1) +
     1          uder(icomp,2,k)*sigma_sl(2) -
     1          uder(icomp,3,k)*sigma_sl(3))
            enddo
         endif
      endif
c
      if (ifdouble.eq.1) then
c
c     See notes (okada_imageanalysis) for explanation of sign
c     flip with j in loop below.
c
         do j = 1,3
         do k = 1,3
            uder(icomp,j,k)=uder(icomp,j,k)*scale
            if (j.eq.3) uder(icomp,j,k) = -uder(icomp,j,k)
         enddo
         enddo
c
         uuu = 0
         do n=1,3
            uuu = uuu+uder(icomp,n,n)
         enddo
c
         do j=1,3
         do k=1,3
            pot = pot + 
     1        rmu*(uder(icomp,j,k)+uder(icomp,k,j))*dipstr(j)*dipvec(k)
         enddo
         enddo
         dotprod = dipstr(1)*dipvec(1) + 
     1             dipstr(2)*dipvec(2) + 
     1             dipstr(3)*dipvec(3) 
         pot = pot + rlam*uuu*dotprod
         if (iffld.eq.1) then
            udertot(4,0,0) = scale*(
     #-3*(2*rr3*dx**6+6*dx**4*rr3**3+4*cd*dx**4*dy**2+6*rr3
     #*dx**4*dy**2+4*rr3**2*dx**4*cd+2*rr3*dx**2*dy**4+4*dx**2*dy
     #**2*rr3**3+3*rr3**4*dx**2*cd+2*dx**2*rr3**5+6*rr3**2*dx**2*
     #cd*dy**2+3*cd*dx**2*dy**4-cd*
     #dy**6-2*rr3**6*cd-4*rr3**2*cd*
     #dy**4-5*rr3**4*cd*dy**2-2*rr3*dy**6-6*dy**4*rr3**
     #3-6*dy**2*rr3**5-2*rr3**7)/(cd+rr3)**4/(dx**2+dy**2+rr3**2)**3)
         udertot(3,1,0) = scale*(
     # 3*(2*cd*dx**4-dx**2*cd*dy**2
     #-4*rr3**3*dx**2-4*rr3*dx**2*dy**2-rr3**2*dx**2*cd-4*rr3*dy**4
     #-8*rr3**3*dy**2-7*rr3**2*cd*dy**2-4*rr3**4*cd
     #-4*rr3**5-3*cd*dy**4)*dx*dy/(cd+rr3)**4/(dx**2+dy**2+rr3**2)**3)
         udertot(3,0,1) = scale*(
     #dx*(2*cd*dx**4-9*rr3*dx**2*dy**2-9*rr3**3*dx**2-4*dx*
     #*2*cd*dy**2-4*rr3**2*dx**2*cd-6*cd*dy**4-9*rr3**5-
     #9*rr3**4*cd-15*rr3**2
     #*cd*dy**2-9*rr3*dy**4-18*rr3**3*dy**2)/
     #(cd+rr3)**3/(dx**2+dy**2+rr3**2)**3)
         udertot(2,2,0) = scale*(
     #-(2*cd*dx**6+2*rr3*dx**6-9*cd
     #*dx**4*dy**2+3*rr3**2*dx**4*cd+2*dx**4*rr3**3-6*rr3*dx**4*dy
     #**2-12*rr3**2*dx**2*cd*dy**2-8*dx**2*dy**2*rr3**3-rr3**4*dx
     #**2*cd-9*cd*dx**2*dy**4-6*rr3*dx**2
     #*dy**4-2*dx**2*rr3**5-2*dy**2*rr3**5-2*rr3**6*cd-rr3**4*cd
     #*dy**2+2*cd*dy**6+2*dy**4*rr3**3+2*rr3
     #*dy**6+3*rr3**2*cd*dy**4-2*rr3**7)/
     #(cd+rr3)**4/(dx**2+dy**2+rr3**2)**3)
         udertot(2,1,1) = scale*(
     #dy*(6*cd*dx**4+6*rr3*dx**4+3*rr3**3*dx**2+6*rr3**2*dx*
     #*2*cd+3*rr3*dx**2*dy**2+4*dx**2*cd*
     #dy**2-3*rr3**4*cd-2*cd*dy**4-5*rr3**
     #2*cd*dy**2-6*rr3**3*dy**2-3*rr3**5-3*rr3*dy**4)/
     #(cd+rr3)**3/(dx**2+dy**2+rr3**2)**3)
         udertot(2,0,2) = scale*(
     #(2*cd*dx**4+4*rr3*dx**4+2*rr3*dx**2*dy**2+3*rr3**2*dx*
     #*2*cd+2*rr3**3*dx**2+dx**2*cd*dy**2
     #-2*rr3**4*cd-3*rr3**2*cd*dy**2-2*rr3*
     #dy**4-4*rr3**3*dy**2-2*rr3**5-cd*dy**4)/
     #(dx**2+dy**2+rr3**2)**3/(cd+rr3)**2)
         udertot(1,3,0) = scale*(
     #-3*(4*rr3*dx**4+3*cd*dx**4+dx**2*cd*
     #dy**2+7*rr3**2*dx**2*cd+4*rr3*dx**2*dy**2+8*rr3**3*dx
     #**2+rr3**2*cd*dy**2-2*cd*dy**4+4*rr3
     #**4*cd+4*rr3**3*dy**2+4*rr3**5)*dx*dy/(cd+rr3)**4/
     #(dx**2+dy**2+rr3**2)**3)
         udertot(1,2,1) = scale*(
     # -(2*cd*dx**4+3*rr3*dx**4+5*rr3**2*dx**2*cd+6*rr3**3*dx**2-
     #4*dx**2*cd*dy**2-3*rr3*dx**2*dy**2
     #+3*rr3**5-6*rr3*dy**4+3*rr3**4*cd-6*cd*
     #dy**4-6*rr3**2*cd*dy**2-3*rr3**3*dy**2)*dx/(cd+rr3)**3/
     #(dx**2+dy**2+rr3**2)**3)
         udertot(1,1,2) = scale*(
     #3*(cd*dx**2+2*rr3*dx**2+2*cd*
     #rr3**2+2*rr3*dy**2+cd*dy**2+2*rr3**3)*dx*dy/(dx**2+dy**2+rr3**
     #2)**3/(cd+rr3)**2)
         udertot(1,0,3) = scale*3/cd**5*dx*rr3
         udertot(0,4,0) = scale*(
     #3*(cd*dx**6+2*rr3*dx**6+6*dx**4*rr3**3+4*rr3**2*dx**4
     #*cd-3*cd*dx**4*dy**2-2*rr3*dx**4*dy*
     #*2-6*rr3**2*dx**2*cd*dy**2+6*dx**2*rr3**5+5*rr3**4*dx**2*s
     #qrt(dx**2+dy**2+rr3**2)-6*rr3*dx**2*dy**4-4*cd*dx**2*dy**4
     #-4*dx**2*dy**2*rr3**3-6*dy**4*rr3**3-2*rr3*dy**6-4*rr3**2*cd
     #*dy**4+2*rr3**6*cd+2*rr3**7-3*rr3**4*cd
     #*dy**2-2*dy**2*rr3**5)/(cd+rr3)**4/(dx**2+dy**2+rr3**2)*
     #*3)
         udertot(0,3,1) = scale*(
     # -(9*rr3*dx**4+6*cd*dx**4+18*rr3**3*dx**2+15*rr3**2*dx
     #**2*cd+9*rr3*dx**2*dy**2+4*dx**2*cd
     #*dy**2+9*rr3**5+9*rr3**4*cd-2*cd*dy*
     #*4+4*rr3**2*cd*dy**2+9*rr3**3*dy**2)*dy/(cd+rr3)**3/
     #(dx**2+dy**2+rr3**2)**3)
         udertot(0,2,2) =  scale*(
     #-(2*rr3*dx**4+cd*dx**4+4*rr3**3*dx**2+3*rr3**2*dx**2*
     #cd-2*rr3*dx**2*dy**2-dx**2*cd*dy**2-
     #2*cd*dy**4-4*rr3*dy**4+2*rr3**5-3*rr3**2*cd
     #*dy**2+2*rr3**4*cd-2*rr3**3*dy**2)/(dx**2+dy**2+rr3**
     #2)**3/(cd+rr3)**2)
         udertot(0,1,3) =  scale*3/cd**5*dy*rr3
         udertot(0,0,4) = -scale*(dx**2+dy**2-2*rr3**2)/cd**5
c
         do l = 1,3
            uuuder(l) = 0
            do n=1,3
               nder(1) = 0
               nder(2) = 0
               nder(3) = 0
               nder(icomp) = nder(icomp)+1
               nder(l) = nder(l)+1
               nder(n) = nder(n)+2
c
c              uder4 should be d_{x_l} (d u_{icomp}^n / d \xi_n)
c
               uder4 = -udertot(nder(1),nder(2),nder(3))
               if (n.eq.3) uder4 = -uder4
               if (l.eq.3) uder4 = -uder4
               uuuder(l) = uuuder(l)+uder4
ccc               write(6,*)' icomp,l,n',icomp,l,n
ccc               write(6,*)' nder',nder
ccc               write(6,*)' uder4 =',uder4
            enddo
         enddo
c
         do l=1,3
         do j=1,3
         do k=1,3
            nder(1) = 0
            nder(2) = 0
            nder(3) = 0
            nder(icomp) = nder(icomp)+1
            nder(l) = nder(l)+1
            nder(j) = nder(j)+1
            nder(k) = nder(k)+1
            uder1 = -udertot(nder(1),nder(2),nder(3))
            uder2 = -udertot(nder(1),nder(2),nder(3))
            if (k.eq.3) uder2 = -uder2
            if (j.eq.3) uder1 = -uder1
            if (l.eq.3) uder1 = -uder1
            if (l.eq.3) uder2 = -uder2
            fld(l) = fld(l) + 
     1               rmu*(uder1+uder2)*dipstr(j)*dipvec(k)
         enddo
         enddo
         fld(l)=fld(l)+rlam*uuuder(l)*dotprod
         enddo
         endif
      endif
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine intker_mindlinc(icomp,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,ztrg,pot,iffld,fld)
c***********************************************************************
c
c     This subroutine computes the displacement due to the 
c     harmonic Mindlin C image solution and all of its
c     derivatives, from which stress and strain can be extracted.
c     (See okada_imageanalysis for details.)
c     It doesn't include the 1/rmu scaling.
c
c     INPUT:
c
c     icomp      desired component of displacement vector
c     rlame(2)   Lame coefficients supplied in the form
c                rlame(1) = rlam, rlame(2) = rmu
c     source(3)  image source location.
c                This source must be the reflected image in the 
c                upper half-space, NOT the original source in the 
c                lower half-space.
c     ifsingle   (the single layer kernel flag)
c                0 means ignore the <<force>> vector 
c                1 means a <<force>> vector is supplied 
c     sigma_sl   force vector
c     ifdouble   (the double layer kernel flag)
c                0 means ignore the <<dislocation>> 
c                1 means a <<dislocation vector>> is supplied 
c                  with surface normal <<dipvec>>
c     dipstr     dislocation vector
c     dipvec     surface orientation
c                
c     ztrg       target location
c     pot        component <<icomp>> of displacement vector
c     fld        -gradient(pot) at ztrg
c-----------------------------------------------------------------------
      implicit none
      integer icomp,ifsingle,ifdouble,iffld,nder(3)
      integer i,j,k,l,n
      real *8 source(3),ztrg(3), sigma_sl(3), delta(3,3)
      real *8 uder(3,3,3),rlame(2),uuuder(3)
      real *8 v(3,3,3)
      real *8 vdertot(0:4,0:4,0:4)
      real *8 dipstr(3),dipvec(3),rr(3)
      real *8 alpha,d3,d5,d7,d9,dotprod,dx,dy,dz,cd
      real *8 rlam,rmu,rr3,scale2c,sgn,uder1,uder2,uder4,uuu
      real *8 v1,v2,v3,v11,v12,v13,v22,v23,v33,scale2
      complex *16 pot,fld(3)
c      
      rlam = rlame(1)
      rmu = rlame(2)
      alpha = (rlam+rmu)/(rlam+2*rmu)
      scale2 = -(2.0d0-alpha)
ccc      write(6,*)' rlam =',rlam
ccc      write(6,*)' rmu =',rmu
ccc      write(6,*)' alpha =',alpha
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
      dx = ztrg(1)-source(1)
      dy = ztrg(2)-source(2)
      dz = ztrg(3)-source(3)
      rr(1) = dx
      rr(2) = dy
      rr(3) = -dz
      rr3 = rr(3)
      cd = sqrt(dx*dx+dy*dy+dz*dz)
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
c
c     Compute three derivatives of harmonic function 
c     V = 1/R
c     with respect to R_1, R_2 and R_3. These are V1,V11,
c     V(1,1,1) etc. We will use these
c     derivatives to compute both the gradient of the SLP
c     and the DLP.
c
c     In this subroutine, it is convenient to first compute
c     uder(i,j,k) = -gradient with respect to R_k of
c     SLP kernel (assuming xi_3 is fixed). 
c     When used to evaluate gradient of SLP below, 
c     uder(i,j,k) needs sign flip for component k=3, since 
c
c     d/dx_1 = d/dR_1, d/dx_2 = d/dR_2, d/dx_3 = - d/dR_3
c
c     It is, however, correctly structured for use in the DLP
c     kernel (if present).
c
      d9 = 1.0d0/(cd**9)
      d7 = 1.0d0/(cd**7)
      d5 = 1.0d0/(cd**5)
      d3 = 1.0d0/(cd**3)
      v1 = -rr(1)*d3
      v2 = -rr(2)*d3
      v3 = -rr(3)*d3
      v11 = -d3 + 3*d5*rr(1)**2
      v12 = 3*d5*rr(1)*rr(2)
      v13 = 3*d5*rr(1)*rr(3)
      v22 = -d3 + 3*d5*rr(2)**2
      v23 = 3*d5*rr(2)*rr(3)
      v33 = -d3 + 3*d5*rr(3)**2
      if ( (ifdouble.eq.1) .or. (iffld.eq.1) ) then
         v(1,1,1) = 9*rr(1)*d5 - 15*d7*rr(1)**3
         v(1,1,2) = 3*rr(2)*d5 - 15*d7*rr(2)*rr(1)**2
         v(1,1,3) = 3*rr(3)*d5 - 15*d7*rr(3)*rr(1)**2
         v(1,2,2) = 3*rr(1)*d5 - 15*d7*rr(1)*rr(2)**2
         v(1,2,3) = - 15*d7*rr(1)*rr(2)*rr(3)
         v(1,3,3) = 3*rr(1)*d5 - 15*d7*rr(1)*rr(3)**2
         uder(1,1,1) = v(1,1,1)
         uder(1,1,2) = v(1,1,2)
         uder(1,2,1) = v(1,1,2)
         v(1,2,1) = v(1,1,2)
         uder(1,1,3) = v(1,1,3)
         uder(1,3,1) = v(1,1,3)
         v(1,3,1) = v(1,1,3)
         uder(1,3,3) = v(1,3,3)
         uder(1,2,3) = v(1,2,3)
         uder(1,3,2) = v(1,2,3)
         v(1,3,2) = v(1,2,3)
         uder(1,2,2) = v(1,2,2)
         if (icomp.ne.1) then
            uder(2,1,1) = v(1,1,2)
            v(2,1,1) = v(1,1,2)
            v(2,2,2)  = 9*rr(2)*d5 - 15*d7*rr(2)**3
            v(2,2,3)  = 3*rr(3)*d5 - 15*d7*rr(3)*rr(2)**2
            v(2,3,3)  = 3*rr(2)*d5 - 15*d7*rr(2)*rr(3)**2
            uder(2,2,2) = v(2,2,2)
            uder(2,1,2) = v(1,2,2)
            v(2,1,2) = v(1,2,2)
            uder(2,2,1) = v(1,2,2)
            v(2,2,1) = v(1,2,2)
            uder(2,2,3) = v(2,2,3)
            uder(2,3,2) = v(2,2,3)
            v(2,3,2) = v(2,2,3)
            uder(2,3,3) = v(2,3,3)
            uder(2,1,3) = v(1,2,3)
            v(2,1,3) = v(1,2,3)
            uder(2,3,1) = v(1,2,3)
            v(2,3,1) = v(1,2,3)
         endif
         if (icomp.eq. 3) then
            uder(3,1,1) = v(1,1,3)
            v(3,1,1) = v(1,1,3)
            uder(3,2,2) = v(2,2,3)
            v(3,2,2) = v(2,2,3)
            v(3,3,3) = 9*rr(3)*d5 - 15*d7*rr(3)**3
            uder(3,3,3) = v(3,3,3)
            uder(3,1,3) = v(1,3,3)
            v(3,1,3) = v(1,3,3)
            uder(3,3,1) = v(1,3,3)
            v(3,3,1) = v(1,3,3)
            uder(3,2,3) = v(2,3,3)
            v(3,2,3) = v(2,3,3)
            uder(3,3,2) = v(2,3,3)
            v(3,3,2) = v(2,3,3)
            uder(3,1,2) = v(1,2,3)
            v(3,1,2) = v(1,2,3)
            uder(3,2,1) = v(1,2,3)
            v(3,2,1) = v(1,2,3)
         endif
         do j = 1,3
            do k = 1,3
               uder(icomp,j,k) = -uder(icomp,j,k)*alpha*source(3)
               if (icomp.eq.3) uder(icomp,j,k) = -uder(icomp,j,k)
            enddo
         enddo
         if (icomp.eq.1) then
            uder(1,3,1)=uder(1,3,1)-scale2*v11
            uder(1,3,2)=uder(1,3,2)-scale2*v12
            uder(1,3,3)=uder(1,3,3)-scale2*v13
         endif
         if (icomp.eq.2) then
            uder(2,3,1)=uder(2,3,1)-scale2*v12
            uder(2,3,2)=uder(2,3,2)-scale2*v22
            uder(2,3,3)=uder(2,3,3)-scale2*v23
         endif
         if (icomp.eq.3) then
            uder(3,1,1)=uder(3,1,1)-scale2*v11
            uder(3,1,2)=uder(3,1,2)-scale2*v12
            uder(3,1,3)=uder(3,1,3)-scale2*v13
            uder(3,2,1)=uder(3,2,1)-scale2*v12
            uder(3,2,2)=uder(3,2,2)-scale2*v22
            uder(3,2,3)=uder(3,2,3)-scale2*v23
         endif
      endif
c
      if (ifsingle.eq.1) then
c
c     compute contribution to SLP for component icomp
c
         scale2c = -(2.0d0-alpha)*sigma_sl(3)
         if (icomp.eq.1) then
            pot = scale2c*v1 + source(3)* 
     1         alpha*(sigma_sl(1)*v11+sigma_sl(2)*v12+sigma_sl(3)*v13)
         else if (icomp.eq.2) then
            pot = scale2c*v2 + source(3)* 
     1         alpha*(sigma_sl(1)*v12+sigma_sl(2)*v22+sigma_sl(3)*v23)
         else 
            pot = scale2*(v1*sigma_sl(1)+v2*sigma_sl(2)) - source(3)* 
     1         alpha*(sigma_sl(1)*v13 + sigma_sl(2)*v23 + 
     1         sigma_sl(3)*v33)
         endif 
c
c
c     compute contribution to -gradient of SLP for component icomp
c
        if (iffld.eq.1) then
           do i=icomp,icomp
           do k=1,3
           sgn = 1.0d0
              if (k.eq.3) sgn = -1.0d0
              do j=1,3
                 fld(k) = fld(k) + 
     1               sgn*uder(i,j,k)*sigma_sl(j)
              enddo
           enddo
           enddo
         endif
      endif
ccc      write(6,*)' after slp'
ccc      write(6,*)(dreal(fld(iii)),iii=1,3)
c
c
      if (ifdouble.eq.1) then
c
c     Compute the double layer kernel.
c
c     Recall that above, we computed
c     uder(i,j,k) = -gradient with respect to R_k of
c     SLP kernel (assuming xi_3 is fixed). 
c     Since
c
c     d/dxi_1 = -d/dR_1, d/dxi_2 = -d/dR_2, d/dxi_3 = -d/dR_3
c
c     uder(i,j,k) can now be interpreted as DLP kernel.  
c
c     However, from C image formula, we clearly need
c     an additional contribution because xi_3
c     appears in definition of kernel and that was held fixed
c     above when we differentiated w.r.t. x_3.
c
         if (icomp.eq.1) then
            uder(1,1,3) = uder(1,1,3) - alpha*v11
            uder(1,2,3) = uder(1,2,3) - alpha*v12
            uder(1,3,3) = uder(1,3,3) - alpha*v13
         else if (icomp.eq.2) then
            uder(2,1,3) = uder(2,1,3) - alpha*v12
            uder(2,2,3) = uder(2,2,3) - alpha*v22
            uder(2,3,3) = uder(2,3,3) - alpha*v23
         else if (icomp.eq.3) then
            uder(3,1,3) = uder(3,1,3) + alpha*v13
            uder(3,2,3) = uder(3,2,3) + alpha*v23
            uder(3,3,3) = uder(3,3,3) + alpha*v33
         endif
c
c     Compute vdertot to hold fourth derivatives of V.
c
      if (iffld.eq.1) then
      vdertot(4,0,0) = 3*(8*dx**4-24*dx**2*dy**2-24*dx**2*rr3**2+
     #3*dy**4+6*dy**2*rr3**2+3*rr3**4)*d9
      vdertot(3,1,0) = 15*dx*dy*(4*dx**2-3*dy**2-3*rr3**2)*d9
      vdertot(2,2,0) = -3*(-27*dx**2*dy**2+4*dx**4+3*dx**2*rr3**2+
     #4*dy**4+3*dy**2*rr3**2-rr3**4)*d9
      vdertot(2,1,1) = 15*dy*rr3*(6*dx**2-dy**2-rr3**2)*d9
      vdertot(3,0,1) = 15*dx*rr3*(4*dx**2-3*dy**2-3*rr3**2)*d9
      vdertot(2,0,2) = -3*(-27*dx**2*rr3**2+4*dx**4+3*dx**2*dy**2+
     #3*dy**2*rr3**2+4*rr3**4-dy**4)*d9
      vdertot(1,3,0) = -15*dx*dy*(-4*dy**2+3*dx**2+3*rr3**2)*d9
      vdertot(1,2,1) = -15*rr3*dx*(-6*dy**2+dx**2+rr3**2)*d9
      vdertot(1,1,2) = -15*dx*dy*(-6*rr3**2+dx**2+dy**2)*d9
      vdertot(1,0,3) = -15*dx*rr3*(-4*rr3**2+3*dx**2+3*dy**2)*d9
      vdertot(0,4,0) = 3*(8*dy**4-24*dx**2*dy**2-24*dy**2*rr3**2+
     #3*dx**4+6*dx**2*rr3**2+3*rr3**4)*d9
      vdertot(0,3,1) = -15*dy*rr3*(-4*dy**2+3*dx**2+3*rr3**2)*d9
      vdertot(0,2,2) = 3*(27*dy**2*rr3**2-3*dx**2*dy**2-4*dy**4-
     #3*dx**2*rr3**2-4*rr3**4+dx**4)*d9
      vdertot(0,1,3) = -15*dy*rr3*(-4*rr3**2+3*dx**2+3*dy**2)*d9
      vdertot(0,0,4) = 3*(8*rr3**4-24*dx**2*rr3**2-24*dy**2*rr3**2+
     #3*dx**4+6*dx**2*dy**2+3*dy**4)*d9
      endif
c
c    compute uuu and uuuder 
c    uuu = sum_n du_{icomp}^n/d\xi_n
c    uuuder(l) = sum_n -d/dx_l du_{icomp}^n/d\xi_n
c
c
         uuu = 0
         do n=1,3
            uuu = uuu+uder(icomp,n,n)
         enddo
c
         if( iffld .eq. 1 ) then
         do l = 1,3
            uuuder(l) = 0
            do n=1,3
               nder(1) = 0
               nder(2) = 0
               nder(3) = 0
               nder(icomp) = nder(icomp)+1
               nder(l) = nder(l)+1
               nder(n) = nder(n)+2
c
c     nder counts number of derivs w.r.t R_1, R_2, R_3.
c
               uder4= -vdertot(nder(1),nder(2),nder(3))*alpha*source(3)
               if (icomp.eq.3) uder4 = -uder4
               uuuder(l) = uuuder(l)+uder4
               if (icomp.eq.3) uuuder(l) = uuuder(l) - scale2*v(n,n,l)
            enddo
            if (icomp.ne.3) uuuder(l) = uuuder(l) - scale2*v(icomp,3,l)
            if (icomp.eq.3) uuuder(l) = uuuder(l) + scale2*v(3,3,l)
            if (icomp.eq.3) uuuder(l) = uuuder(l) + alpha*v(3,3,l)
            if (icomp.ne.3) uuuder(l) = uuuder(l) - alpha*v(icomp,3,l)
            if (l.ne.3) uuuder(l) = -uuuder(l)
         enddo
         endif
c
c     Compute contribution to displacement component.
c
         dotprod = dipstr(1)*dipvec(1) + 
     1             dipstr(2)*dipvec(2) + 
     1             dipstr(3)*dipvec(3) 
         do j=1,3
         do k=1,3
            pot = pot + 
     1      rmu*(uder(icomp,j,k)+uder(icomp,k,j))*dipstr(j)*dipvec(k)
         enddo
         enddo
         pot = pot + rlam*uuu*dotprod
        if (iffld.eq.1) then
c
c     Compte contributions to -gradient of displacement compoment.
c
         do l=1,3
         do j=1,3
         do k=1,3
            nder(1) = 0
            nder(2) = 0
            nder(3) = 0
            nder(icomp) = nder(icomp)+1
            nder(l) = nder(l)+1
            nder(j) = nder(j)+1
            nder(k) = nder(k)+1
c
c     nder counts number of derivs w.r.t R_1, R_2, R_3.
c     uder1 will hold d/dR_l du_{icomp}^j/d\xi_k
c     uder2 will hold d/dR_l du_{icomp}^k/d\xi_j
c
c     lots of if statements to account for the numerous
c     special case contributions to kernel.
c
            uder1 = -vdertot(nder(1),nder(2),nder(3))*alpha*source(3)
            uder2 = -vdertot(nder(1),nder(2),nder(3))*alpha*source(3)
            if (icomp.eq.3) uder1 = -uder1
            if (icomp.eq.3) uder2 = -uder2
            if ((icomp.ne.3).and.(j.eq.3)) uder1 = uder1 -
     1                         scale2*v(icomp,k,l)
            if ((icomp.ne.3).and.(k.eq.3)) uder2 = uder2 -
     1                         scale2*v(icomp,j,l)
            if ((icomp.eq.3).and.(j.ne.3)) 
     1                     uder1 = uder1 - scale2*v(j,k,l)
            if ((icomp.eq.3).and.(k.ne.3)) 
     1                     uder2 = uder2 - scale2*v(j,k,l)
c
            if ((icomp.ne.3).and.(k.eq.3)) uder1 = uder1 -
     1                         alpha*v(icomp,j,l)
            if ((icomp.ne.3).and.(j.eq.3)) uder2 = uder2 -
     1                         alpha*v(icomp,k,l)
            if ((icomp.eq.3).and.(k.eq.3)) uder1 = uder1 +
     1                         alpha*v(icomp,j,l)
            if ((icomp.eq.3).and.(j.eq.3)) uder2 = uder2 +
     1                         alpha*v(icomp,k,l)
c
c     want -d/dx_l, so flip sign for first and second
c     component. (d/dR_3 is already - d/dx_3.)
c
            if (l.ne.3) uder1 = -uder1
            if (l.ne.3) uder2 = -uder2
            fld(l) = fld(l) + 
     1               rmu*(uder1+uder2)*dipstr(j)*dipvec(k)
           fld(l)=fld(l)+rlam*delta(j,k)*uuuder(l)*dipstr(j)*dipvec(k)
         enddo
         enddo
         enddo
        endif
      endif
ccc      write(6,*)' after dlp'
ccc      write(6,*)(dreal(fld(iii)),iii=1,3)
      return
      end

