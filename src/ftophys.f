C***********************************************************************
      subroutine ftophys(mexpf,nexptot,nlambs,rlams,numfour,numphys,
     1                      nthmax,mexpphys,nexptotp,fexpe,nexpe)
C***********************************************************************
C
C     This subroutine converts the Fourier expansion of the
C     exponential moment function into a physical plane wave expansion.
C     nodes.
C
C     INPUT:
C
C     mexpf     =  Fourier coefficients of the function 
C                  MEXP(lambda,alpha) for discrete lambda values. 
C                  They are ordered as follows:
C
C                 mexpf(1,...,numfour(1)) = Fourier modes for lambda_1
C                 mexpf(numfour(1)+1,...,numfour(2)) = Fourier modes
C                                              for lambda_2
C                 etc.
C
C     nexptot   =  length of mexpf
C     nlambs    =  number of discretization pts. in lambda integral
C     rlams     =  discretization points in lambda integral.
C     numfour   =  number of Fourier modes in the expansion
C                      of the function MEXP(\lambda_j,\alpha)
C     numphys   =  number of desired samples in \alpha variable 
C                      of the function MEXP(\lambda_j,\alpha)
C     nthmax    =  max_j numfour(j)
C     nexptotp  =  length of output array mexpphys
C     fexpe     =  precomputed array of exponentials needed for
C                   Fourier series evaluation
C     nexpe     =  length of fexpe
C
C     OUTPUT:
C
C     mexpphys  =  Discrete values of the moment function 
C                   MEXP(\lambda,\alpha), ordered as follows.
C
C         mexpphys(1),...,mexpphys(numphys(1)) = MEXP(\lambda_1,0),..., 
C              MEXP(\lambda_1, 2*PI*(numphys(1)-1)/numphys(1)).
C         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
C              MEXP(\lambda_2,0),...,
C                  MEXP(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
C         etc.
C
C------------------------------------------------------------
      implicit none
      integer nlambs,numfour(nlambs),numphys(nlambs),nthmax
      integer nexpe,nexptot,nexptotp,mm,nftot,nptot,nexte,i,ival
      complex *16 mexpf(nexptot)
      complex *16 mexpphys(nexptotp)
      complex *16 fexpe(nexpe)
      complex *16 ctmp
      real *8     rlams(nlambs)
      real *8     sgn
C
      nftot = 0
      nptot  = 0
c
      nexte = 1
      do i=1,nlambs
        do ival=1,numphys(i)
           mexpphys(nptot+ival) = mexpf(nftot+1)
           sgn = -1
           do mm = 2,numfour(i)
              ctmp = fexpe(nexte)*mexpf(nftot+mm)
              mexpphys(nptot+ival) = mexpphys(nptot+ival) +
     1                ctmp + dconjg(ctmp)*sgn
              nexte = nexte + 1
              sgn = -sgn
           enddo
        enddo
        nftot = nftot+numfour(i)
        nptot = nptot+numphys(i)
      enddo
      return
      end
c
