C***********************************************************************
      subroutine phystof(mexpf,nexptot,nlambs,rlams,numfour,numphys,
     1                      nthmax,mexpphys,nexptotp,fexpback,nexpback)
c***********************************************************************
C
C     This subroutine converts the discretized exponential moment function
C     into its Fourier expansion.
C
C     On OUTPUT:
C
C     mexpf(*):     Fourier coefficients of the function 
C                   M(lambda,alpha) for discrete lambda values. 
C                   They are ordered as follows:
C
C               mexpf(1,...,numfour(1)) = Fourier modes for lambda_1
C               mexpf(numfour(1)+1,...,numfour(2)) = Fourier modes
C                                              for lambda_2
C               etc.
C     INPUT:
C
C     nexptot:       length of mexpf
C     nlambs:        number of discretization pts. in lambda integral
C     rlams(nlambs): discretization points in lambda integral.
C     numfour(j):    number of Fourier modes in the expansion
C                      of the function M(\lambda_j,\alpha)
C     numphys(j):    number of samples of the function
C                      M(\lambda_j,\alpha) in \alpha variable
C     nthmax =      max_j numfour(j)
C     mexpphys(*)   Discrete values of the moment function 
C                   M(\lambda,\alpha), ordered as follows.
C
C            mexpphys(1),...,mexpphys(numphys(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
C            mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
C            etc.
C
C     nexptotp =    length of mexpphys
C     fexpback =    precomputed array of exponentials needed for
C                  Fourier series evaluation
C     nexpback =    length of fexpback
C
C
C------------------------------------------------------------
      implicit none
      integer     nlambs,numfour(nlambs),numphys(nlambs),nthmax
      integer     nexptot,nexptotp,nexpback,nftot,nptot
      integer     next,i,mm,nalpha,ival
      real *8     rlams(nlambs)
      real *8     halpha,done,pi 
      complex *16 mexpf(nexptot)
      complex *16 mexpphys(nexptotp)
      complex *16 fexpback(nexpback)
C------------------------------------------------------------
C
      done=1.0D0
      pi=datan(done)*4
C
      nftot = 0
      nptot  = 0
      next  = 1
      do i=1,nlambs
  	 nalpha = numphys(i)
         halpha=2*pi/nalpha
         mexpf(nftot+1) = 0.0d0
         do ival=1,nalpha
            mexpf(nftot+1) = mexpf(nftot+1) + mexpphys(nptot+ival) 
         enddo
         mexpf(nftot+1) = mexpf(nftot+1)/nalpha
         do mm = 2,numfour(i)
            mexpf(nftot+mm) = 0.0d0
            do ival=1,nalpha
              mexpf(nftot+mm) = mexpf(nftot+mm) +
     1          fexpback(next)*mexpphys(nptot+ival)
                next = next+1
            enddo
            mexpf(nftot+mm) = mexpf(nftot+mm)/nalpha
         enddo
         nftot = nftot+numfour(i)
         nptot = nptot+numphys(i)
      enddo
      return
      end
c
