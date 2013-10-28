C***********************************************************************
      subroutine mkfexp(nlambs,numfour,numphys,fexpe,nexpe,
     1           fexpback,nexpback)
C***********************************************************************
C
C     This subroutine computes the tables of exponentials needed
C     for mapping from Fourier to physical domain and back.
C     In order to minimize storage, they are organized in a 
C     one-dimenional array corresponding to the order in which they
C     are accessed by subroutines FTOPHYS and PHYSTOF.
C    
c-----------------------------------------------------------------------
C     INPUT:
C   
C     nlambs    =  number of discretization pts in lambda integral
C     numfour   =  number of Fourier modes in the expansion
C                      of the function MEXP(\lambda_j,\alpha)
C     numphys   =  number of desired samples in \alpha variable 
C                      of the function MEXP(\lambda_j,\alpha)
C     nexpe     =  length of fexpe
C     nexpback  =  length of fexpback
C
C     nexpe, nexpback should be computed in prior call to 
C     subroutine getfexplengths.
C
C     OUTPUT:
C
C     fexpe     =  array of exponentials needed to map from
C                   Fourier representation to physical expansion
C     fexpback  =  array of exponentials needed to map from
C                   physical representation to Fourier expansion
C
c-----------------------------------------------------------------------
      implicit none
      integer nlambs,numphys(nlambs),numfour(nlambs)
      integer nexpe,nexpback
      integer nexte,next,i,nalpha,j,mm
      real *8 halpha,alpha,pi
      complex *16 fexpe(nexpe)
      complex *16 fexpback(nexpback)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
C
      pi = 4*datan(1.0d0)
      nexte = 1
      do i=1,nlambs
	 nalpha = numphys(i)
         halpha=2*pi/nalpha
         do j=1,nalpha
            alpha=(j-1)*halpha
	    do mm = 2,numfour(i)
               fexpe(nexte)  = cdexp(ima*(mm-1)*alpha)
	       nexte = nexte + 1
            enddo
         enddo
      enddo
c
      next = 1
      do i=1,nlambs
	 nalpha = numphys(i)
         halpha=2*pi/nalpha
	 do mm = 2,numfour(i)
            do j=1,nalpha
               alpha=(j-1)*halpha
               fexpback(next)  = cdexp(-ima*(mm-1)*alpha)
	       next = next + 1
            enddo
         enddo
      enddo
      return
      end
C
C
C
C***********************************************************************
      subroutine getfexplengths(nlambs,numfour,numphys,nexpe,nexpback)
C***********************************************************************
C     This subroutine computes the lengths of the 
C     tables of exponentials needed
C     for mapping from Fourier to physical domain and back.
C     In order to minimize storage, they are organized in a 
C     one-dimenional array corresponding to the order in which they
C     are accessed by subroutines FTOPHYS and PHYSTOF.
C    
c-----------------------------------------------------------------------
C     INPUT:
C   
C     nlambs    =  number of discretization pts in lambda integral
C     numfour   =  number of Fourier modes in the expansion
C                      of the function MEXP(\lambda_j,\alpha)
C     numphys   =  number of desired samples in \alpha variable 
C                      of the function MEXP(\lambda_j,\alpha)
C     OUTPUT:
C
C     nexpe     =  length of fexpe
C     nexpback  =  length of fexpback
C-----------------------------------------------------------------------
      implicit none
      integer nlambs,numphys(nlambs),numfour(nlambs)
      integer nexpe,nexpback,nalpha,i,j,mm
c
      nexpe = 1
      do i=1,nlambs
	 nalpha = numphys(i)
         do j=1,nalpha
	    do mm = 2,numfour(i)
	       nexpe = nexpe + 1
            enddo
         enddo
      enddo
c
      nexpback = 1
      do i=1,nlambs
	 nalpha = numphys(i)
	 do mm = 2,numfour(i)
            do j=1,nalpha
	       nexpback = nexpback + 1
            enddo
         enddo
      enddo
      return
      end
c
c
c
      subroutine mkexps(rlams,nlambs,numphys,nexptotp,xs,ys,zs)
      implicit none
      integer   nlambs,numphys(nlambs)
      integer   nexptotp,nl,mth,ntot,ncurrent
      real *8     rlams(nlambs),u,hu,pi
      complex *16 ima
      complex *16 xs(3,nexptotp)
      complex *16 ys(3,nexptotp)
      real *8 zs(3,nexptotp)
      data ima/(0.0d0,1.0d0)/
c
C     This subroutine computes the tables of exponentials needed
C     for translating plane wave representations of harmonic
C     functions.
C
C     U   = \int_0^\infty e^{-\lambda z}
C           \int_0^{2\pi} e^{i\lambda(x cos(u)+y sin(u))}
C           MEXPPHYS(lambda,u) du dlambda
C
C     MEXPPHYS(*):  Discrete values of the moment function 
C                   M(\lambda,u), ordered as follows.
C
C         MEXPPHYS(1),...,MEXPPHYS(NUMPHYS(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*PI*(NUMPHYS(1)-1)/NUMPHYS(1)).
C         MEXPPHYS(NUMPHYS(1)+1),...,MEXPPHYS(NUMPHYS(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*PI*(NUMPHYS(2)-1)/NUMPHYS(2)).
C         etc.
C
C     INPUT:
C
C     rlams(nlambs)  discretization points in lambda integral 
C     nlambs         number of discret. pts. in lambda integral
C     numphys(j)     number of nodes in u integral needed 
C                    for corresponding lambda =  lambda_j. 
C     nexptotp       sum_j NUMPHYS(j)
C
C     OUTPUT:
C
C     xs(1,nexptotp)   e^{i \lambda_j (cos(u_k)}  in above ordering
C     xs(2,nexptotp)   e^{i \lambda_j (2 cos(u_k)}  in above ordering.
C     xs(3,nexptotp)   e^{i \lambda_j (3 cos(u_k)}  in above ordering.
C     ys(1,nexptotp)   e^{i \lambda_j (sin(u_k)}  in above ordering.
C     ys(2,nexptotp)   e^{i \lambda_j (2 sin(u_k)}  in above ordering.
C     ys(3,nexptotp)   e^{i \lambda_j (3 sin(u_k)}  in above ordering.
C     zs(1,nexptotp)   e^{-\lambda_j}     in above ordering.
C     zs(2,nexptotp)   e^{-2 \lambda_j}   in above ordering. 
C     zs(3,nexptotp)   e^{-3 \lambda_j}   in above ordering. 
C------------------------------------------------------------
C      
C     Loop over each lambda value 
C
      pi = 4*datan(1.0d0)
      ntot = 0
      do nl = 1,nlambs
         hu=2*pi/numphys(nl)
         do mth = 1,numphys(nl)
            u = (mth-1)*hu
            ncurrent = ntot+mth
            zs(1,ncurrent) = dexp( -rlams(nl) )
            zs(2,ncurrent) = dexp( - 2.0d0*rlams(nl) )
            zs(3,ncurrent) = dexp( - 3.0d0*rlams(nl) )
            xs(1,ncurrent) = cdexp(ima*rlams(nl)*dcos(u))
            xs(2,ncurrent) = cdexp(ima*rlams(nl)*2.0d0*dcos(u))
            xs(3,ncurrent) = cdexp(ima*rlams(nl)*3.0d0*dcos(u))
            ys(1,ncurrent) = cdexp(ima*rlams(nl)*dsin(u))
            ys(2,ncurrent) = cdexp(ima*rlams(nl)*2.0d0*dsin(u))
            ys(3,ncurrent) = cdexp(ima*rlams(nl)*3.0d0*dsin(u))
         enddo
         ntot = ntot+numphys(nl)
      enddo
      return
      end

