C***********************************************************************
      subroutine xyzshift(mexpphys,rlams,nlams,numphys,mexppnew,
     1                  nexptotp,xdis,ydis,zdis)    
C***********************************************************************
c
C     This subroutine shifts a (downward) exponential expansion in
C     the x and y and z directions (XDIS,YDIS,ZDIS).
C
C     U   = \int_0^\infty e^{\lambda z}
C           \int_0^{2\pi} e^{i\lambda(x cos(u)+y sin(u))}
C           M(lambda,u) dalpha du
C
C     INPUT:
C
C     mexpphys =  discrete values of the moment function 
C                   M(\lambda,u), ordered as follows.
C
C         mexpphys(1),...,mexpphys(numphys(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
C         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
C                     etc.
C     rlams    =  discretization points in lambda integral 
C     nlams    =  number of discretization pts. in lambda integral
C     numphys  =  array of number of points needed in expansion
C                 of alpha variable for lambda_j. 
C     nexptotp =  size of mexpphys, mexpphysnew
C     xdis:         distance to be shifted in x direction.
C     ydis:         distance to be shifted in y direction.
C     zdis:         distance to be shifted in z direction.
C
C     OUTPUT:
C
C     mexppnew =  discrete values of the shifted moment function 
C                  ordered as for mexpphys.
C
C------------------------------------------------------------
      implicit none
      integer  nlams,nexptotp,ntot,nl,mth,ncurrent
      integer  numphys(nlams)
      real *8  rlams(nlams)
      complex *16 mexpphys(nexptotp)
      complex *16 mexppnew(nexptotp)
      complex *16 rmul,ima
      real *8  hu,pi,xdis,ydis,zdis,u,rmul0
      data ima/(0.0d0,1.0d0)/
C-------------------------------------------------------------------
C      
C     Loop over each lambda value 
C
      pi = 4*datan(1.0d0)
      ntot = 0
      do nl = 1,nlams
         hu=2*pi/numphys(nl)
         rmul0 = dexp(rlams(nl)*zdis)
         do mth = 1,numphys(nl)
            u = (mth-1)*hu
            rmul = rmul0*
     1      cdexp(ima*rlams(nl)*(xdis*dcos(u)+ydis*dsin(u)))
            ncurrent = ntot+mth
            mexppnew(ncurrent) = mexpphys(ncurrent)*rmul
         enddo
         ntot = ntot+numphys(nl)
      enddo
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine xyzshift_add(mexpphys,rlams,nlams,numphys,mexppnew,
     1                  nexptotp,xdis,ydis,zdis)    
C***********************************************************************
c
C     This subroutine shifts a (downward) exponential expansion in
C     the x and y and z directions (XDIS,YDIS,ZDIS).
C
C     U   = \int_0^\infty e^{\lambda z}
C           \int_0^{2\pi} e^{i\lambda(x cos(u)+y sin(u))}
C           M(lambda,u) dalpha du
C
C     INPUT:
C
C     mexpphys =  discrete values of the moment function 
C                   M(\lambda,u), ordered as follows.
C
C         mexpphys(1),...,mexpphys(numphys(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
C         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
C                     etc.
C     rlams    =  discretization points in lambda integral 
C     nlams    =  number of discretization pts. in lambda integral
C     numphys  =  array of number of points needed in expansion
C                 of alpha variable for lambda_j. 
C     nexptotp =  size of mexpphys, mexpphysnew
C     xdis:         distance to be shifted in x direction.
C     ydis:         distance to be shifted in y direction.
C     zdis:         distance to be shifted in z direction.
C
C     OUTPUT:
C
C     mexppnew =  discrete values of the shifted moment function 
C                  ordered as for mexpphys.
C
C------------------------------------------------------------
      implicit none
      integer  nlams,nexptotp,ntot,nl,mth,ncurrent
      integer  numphys(nlams)
      real *8  rlams(nlams)
      complex *16 mexpphys(nexptotp)
      complex *16 mexppnew(nexptotp)
      complex *16 rmul,ima
      real *8  hu,pi,xdis,ydis,zdis,u,rmul0
      data ima/(0.0d0,1.0d0)/
C------------------------------------------------------------
C      
C     Loop over each lambda value 
C
      pi = 4*datan(1.0d0)
      ntot = 0
      do nl = 1,nlams
         hu=2*pi/numphys(nl)
         rmul0 = dexp(rlams(nl)*zdis)
         do mth = 1,numphys(nl)
            u = (mth-1)*hu
            rmul = rmul0*
     1      cdexp(ima*rlams(nl)*(xdis*dcos(u)+ydis*dsin(u)))
            ncurrent = ntot+mth
            mexppnew(ncurrent) = mexppnew(ncurrent)+
     1          mexpphys(ncurrent)*rmul
         enddo
         ntot = ntot+numphys(nl)
      enddo
      return
      end
C
C***********************************************************************
      subroutine xyzshift_add_fast(mexpphys,nexptotp,mexppnew,
     1           netshift,xshift,yshift,zshift,ix,iy,iz)    
C***********************************************************************
c
C     This subroutine shifts a (downward) exponential expansion in
C     the x and y and z directions (XDIS,YDIS,ZDIS).
C     This is an accelerated routine, using precomputed translation
C     operators in teh x,y,z directions.
C
C     U   = \int_0^\infty e^{\lambda z}
C           \int_0^{2\pi} e^{i\lambda(x cos(u)+y sin(u))}
C           M(lambda,u) dalpha du
C
C     INPUT:
C
C     mexpphys =  discrete values of the moment function 
C                   M(\lambda,u), ordered as follows.
C
C         mexpphys(1),...,mexpphys(numphys(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
C         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
C                     etc.
C     nexptotp    size of mexpphys, mexpphysnew
C     netshift    comlpex work array 
C     xshift:     precomputed weights for xshift
C     yshift:     precomputed weights for yshift
C     zshift:     precomputed weights for zshift
c     ix          integer number of boxes translating in x direction
c     iy          integer number of boxes translating in y direction
c     iz          integer number of boxes translating in z direction
C
C     OUTPUT:
C
C     mexppnew =  discrete values of the shifted moment function 
C                  ordered as for mexpphys.
C
C------------------------------------------------------------
      implicit none
      integer  nexptotp,ix,iy,iz,i
      real *8 zshift(3,nexptotp)
      complex *16 mexpphys(nexptotp)
      complex *16 mexppnew(nexptotp)
      complex *16 netshift(nexptotp)
      complex *16 xshift(3,nexptotp)
      complex *16 yshift(3,nexptotp)
C--------------------------------------------------------------
C      
C     Build translation operator from xshift,yshift,zshift
C     according to ix,iy,iz. 
C
      if (ix.eq.0) then
         do i = 1,nexptotp
            netshift(i) = zshift(iz,i)
         enddo
      else if (ix.gt.0) then
         do i = 1,nexptotp
            netshift(i) = zshift(iz,i)*xshift(ix,i)
         enddo
      else
         do i = 1,nexptotp
            netshift(i) = zshift(iz,i)/xshift(abs(ix),i)
         enddo
      endif
c
      if (iy.gt.0) then
         do i = 1,nexptotp
            netshift(i) = netshift(i)*yshift(iy,i)
         enddo
      else if (iy.lt.0) then
         do i = 1,nexptotp
            netshift(i) = netshift(i)/yshift(abs(iy),i)
         enddo
      endif
C      
C     Loop over each plane wave coefficient
C
      do i = 1,nexptotp
            mexppnew(i) = mexppnew(i)+
     1          mexpphys(i)*netshift(i)
      enddo
      return
      end
