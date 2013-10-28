C***********************************************************************
      subroutine mpoletoexp(mpole,nterms,nlambs,numtets,
     1                      nexptot,mexpdown,rlsc)
c***********************************************************************
C
C     This subroutine converts a multipole expansion mpole into the
C     corresponding exponential moment function mexpdown.
C
C     U(x,y,z) = \sum_{n=0}^{nterms} \sum_{m=-n,n}
C                mpole(n,m) Y_n^m(cos theta) e^{i m \phi}/r^{n+1}
C
C              = (1/2pi) \int_0^\infty e^{\lambda z}
C                \int_0^{2\pi} e^{-i\lambda(xcos(alpha)+ysin(alpha))}
C                mexpdown(lambda,alpha) dalpha dlambda
C
C     for -z direction.
C
C-----------------------------------------------------------------------
C     NOTE: The multipole expansion is assumed to have been rescaled
C           so that the box containing sources has unit dimension.
C
C     NOTE: We only use MPOLE(n,m) for n,m >= 0, since MPOLE(n,-m)=  
C           DCONJG(MPOLE(n,m)) for real-valued output.
C           Since we store the exponential
C           moment function in the Fourier domain (w.r.t. the alpha
C           variable), we compute
C
C       M_lambda(m) = (i)**m \sum_{n=m}^N c(n,m) MPOLE(n,m) lambda^n
C
C           for m >= 0 only, where c(n,m) = 1/sqrt((n+m)!(n-m)!).
C
C       For possible future reference, it should be noted that
C       it is NOT TRUE that m_lamb(-m) = dconjg(m_lamb(m)).
C       Inspection of the integral formula for Y_n^{-m} shows that
C       m_lamb(-m) = dconjg(m_lamb(m)) * (-1)**m.
C-----------------------------------------------------------------------
C
C     INPUT:
C
C     mpole:    the multipole expansion
C     nterms:    order of the multipole expansion
C     nlambs:    number of discretization pts in lambda integral
C     numtets:   number of Fourier modes needed in expansion
C                    of alpha variable for each lambda value.
C     nexptot:   sum_j numtets(j)
C     rlsc   :   precomputed array of coefficients
C
C
C     OUTPUT:
C
C     mexpdown: Fourier coefficients of the function 
C                     mexp(lambda,alpha) for successive discrete 
C                     lambda values. They are ordered as follows:
C
C                 mexpf(1,...,numtets(1)) = Fourier modes
C                             for lambda_1
C                 mexpf(numtets(1)+1,...,numtets(2)) = Fourier modes
C                             for lambda_2
C                 etc.
C------------------------------------------------------------
      implicit none
      integer     nterms,nl,mth,ncurrent,nm,ntot
      integer     nlambs,numtets(nlambs),nexptot
      real *8     sgn
      real *8     rlsc(nlambs,0:nterms,0:nterms)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mexpdown(nexptot)
      complex *16 zeyep,ztmp1,ztmp2
C------------------------------------------------------------
C      
C     Loop over multipole order to generate mexpdown values.
C
      ntot = 0
      do nl = 1,nlambs
         sgn = -1.0D0
         zeyep = 1.0d0
         do mth = 0,numtets(nl)-1
            ncurrent = ntot+mth+1
            ztmp1 = 0.0d0
            ztmp2 = 0.0d0
            sgn = -sgn
            do nm = mth,nterms,2
               ztmp1 = ztmp1 + 
     1         rlsc(nl,nm,mth)*mpole(nm,mth)
            enddo
            do nm = mth+1,nterms,2
               ztmp2 = ztmp2 + 
     1         rlsc(nl,nm,mth)*mpole(nm,mth)
            enddo
            mexpdown(ncurrent) = (ztmp1 + ztmp2)*zeyep
            zeyep = -zeyep*dcmplx(0.0d0,1.0d0)
         enddo
         ntot = ntot+numtets(nl)
      enddo
      return
      end
c
C***********************************************************************
      subroutine rlscini(rlsc,nlambs,rlams,nterms)
C***********************************************************************
C     See subroutine mpoletoexp for explanation.
C
C     INPUT:
C
C     nlambs: number of discretization pts in lambda integral
C     rlams:  discretization pts in lambda integral
C     nterms: order of multipole expansion
C
C
C     OUTPUT:
C
C     rlsc:   mapping coefficients 
C------------------------------------------------------------
      implicit none
      integer  nl,j,k,nlambs,nterms
      real *8  rmul
      real *8  rlsc(nlambs,0:nterms,0:nterms)
      real *8  rlams(nlambs)
      real *8, allocatable :: rlampow(:)
      real *8, allocatable :: facts(:)
c
      allocate(facts(0:2*nterms))
      allocate(rlampow(0:nterms))
c
      facts(0) = 1.0d0
      do j = 1,2*nterms
	 facts(j) = facts(j-1)*dsqrt(j+0.0d0)
      enddo
c
      do nl = 1,nlambs
c
c     Compute powers of lambda(nl)
c
         rlampow(0) = 1.0D0
         rmul = -rlams(nl)
         do j = 1,nterms
            rlampow(j) = rlampow(j-1)*rmul
         enddo
         do j = 0,nterms
            do k = 0,j
               rlsc(nl,j,k) = rlampow(j)/(facts(j-k)*facts(j+k))
            enddo
         enddo
      enddo
      return
      end
C
