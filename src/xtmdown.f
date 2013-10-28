C***********************************************************************
      subroutine exptolocald(local,nterms,rlams,whts,nlambs,numtets,
     1                      nthmax,nexptot,mexpdown,boxsize,scale,
     1                      rlampow,facts,zeye)
c***********************************************************************
C
C     This subroutine converts the Fourier representation of two
C     exponential moment functions into a local multipole expansion 
C     (with respect to the same box center). 
C
C     U(x,y,z) = \int_0^\infty e^{\lambda z}
C                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
C                mexpdown(lambda,alpha) dalpha dlambda
C
C              = \sum_{n=0}^{nterms} \sum_{m=-n,n}
C                LOCAL(n,m) Y_n^m(cos theta) e^{i m \phi} r^{n}
C
C     INPUT:
C
C     nterms      order of local expansion
C     rlams       discretization points in lambda integral 
C     whts        quadrature weights in lambda integral 
C     nlambs      number of discretization pts. in lambda integral 
C     numtets(j)  number of Fourier modes in expansion of alpha
C                 variable for lambda_j.
C     nthmax      max_j numtets(j)
C     nexptot     sum_j numtets(j)
C     mexpdown    Fourier coefficients of the function
C                     MEXPDOWN for discrete lambda
C                     values. They are ordered as follows:
C
C                 mexpdown(1,...,numtets(1)) = Fourier modes
C                             for lambda_1
C                 mexpdown(numtets(1)+1,...,numtets(2)) = Fourier modes
C                             for lambda_2
C                 etc.
C
C     boxsize     dimension of box (recall that quadrature rule
C                 requires unit box dimensions so length scales need
C                 to be recalibrated).
C     scale       desired scale parameter for local expansion
C     rlampow     real work array
C     facts       real work array
C     zeye        complex work array
C
C     OUTPUT:
C
C     local(0:nterms,-nterms:nterms): output local expansion 
C
C------------------------------------------------------------
      implicit none
      integer     nlambs,nterms,nexptot,nthmax
      integer     i,nm,mth,ntot,nl,mmax,ncurrent
      integer     numtets(nlambs)
      real *8     rlams(nlambs)
      real *8     whts(nlambs)
      real *8     rlampow(0:nterms)
      real *8     facts(0:2*nterms)
      real *8     boxsize,rmul,scale,rscale
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 mexpdown(nexptot)
      complex *16 zeye(0:nthmax)
c-----------------------------------------------------------
      facts(0) = 1.0d0
      do i = 1,2*nterms
	 facts(i) = facts(i-1)*dsqrt(i+0.0d0)
      enddo
c
C     compute necessary powers of -i
C
      zeye(0) = 1.0D0
      do i = 1,nthmax
         zeye(i) = zeye(i-1)*dcmplx(0.0d0,-1.0d0)
      enddo
c
C     initialize local expansion
C
      do nm = 0,nterms
         do mth = -nterms,nterms
            local(nm,mth) = dcmplx(0.0d0,0.0d0)
         enddo
      enddo
C
C     Loop over lambdas computing contributions to 
C     local expansion
C
      ntot = 1
      do nl = 1,nlambs
C
C     compute powers of lambda_nl
C
         rlampow(0) = whts(nl)
         rmul = rlams(nl)/(boxsize*scale)
         do nm = 1,nterms
            rlampow(nm) = rlampow(nm-1)*rmul
         enddo
c
C     add contributions to LOCAL expansion.
c
         do nm = 0,nterms
            mmax = numtets(nl)-1
            if (mmax.gt.nm) mmax = nm
	    rmul = rlampow(nm)
            do mth = 0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth) + rmul*
     1          mexpdown(ncurrent)
            enddo
         enddo
         ntot = ntot+numtets(nl)
      enddo
c
C     scale the expansions according to formula
C
      do nm = 0,nterms
         do mth = 0,nm
            rscale = 1.0d0/(facts(nm-mth)*facts(nm+mth)) 
            local(nm,mth) = local(nm,mth)*zeye(mth)*rscale/boxsize
         enddo
         do mth = 1,nm
            local(nm,-mth) = dconjg(local(nm,mth))
         enddo
      enddo
      return
      end
C
C
C
C***********************************************************************
      subroutine exptolocald_add(local,nterms,rlams,whts,nlambs,numtets,
     1                      nthmax,nexptot,mexpdown,localtmp,boxsize,
     1                      scale,rlampow,facts,zeye)
c***********************************************************************
C     This subroutine converts the Fourier representation of two
C     exponential moment functions into a local multipole expansion 
C     (with respect to the same box center). 
C
C     U(x,y,z) = \int_0^\infty e^{\lambda z}
C                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
C                mexpdown(lambda,alpha) dalpha dlambda
C
C              = \sum_{n=0}^{nterms} \sum_{m=-n,n}
C                LOCAL(n,m) Y_n^m(cos theta) e^{i m \phi} r^{n}
C
C     INPUT:
C
C     nterms      order of local expansion
C     rlams       discretization points in lambda integral 
C     whts        quadrature weights in lambda integral 
C     nlambs      number of discretization pts. in lambda integral 
C     numtets(j)  number of Fourier modes in expansion of alpha
C                 variable for lambda_j.
C     nthmax      max_j numtets(j)
C     nexptot     sum_j numtets(j)
C     mexpdown    Fourier coefficients of the function
C                     MEXPDOWN for discrete lambda
C                     values. They are ordered as follows:
C
C                 mexpdown(1,...,numtets(1)) = Fourier modes
C                             for lambda_1
C                 mexpdown(numtets(1)+1,...,numtets(2)) = Fourier modes
C                             for lambda_2
C                 etc.
C
C     boxsize     dimension of box (recall that quadrature rule
C                 requires unit box dimensions so length scales need
C                 to be recalibrated).
C     scale       desired scale parameter for local expansion
C     rlampow     real work array
C     facts       real work array
C     zeye        complex work array
C
C     OUTPUT:
C
C     local(0:nterms,-nterms:nterms): output local expansion 
C
C------------------------------------------------------------
      implicit none
      integer     nlambs,nterms,nexptot,nthmax
      integer     i,nm,mth,ntot,nl,mmax,ncurrent
      integer     numtets(nlambs)
      real *8     rlams(nlambs)
      real *8     whts(nlambs)
      real *8     rlampow(0:nterms)
      real *8     facts(0:2*nterms)
      real *8     boxsize,rmul,scale,rscale
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 localtmp(0:nterms,-nterms:nterms)
      complex *16 mexpdown(nexptot)
      complex *16 zeye(0:nthmax)
C------------------------------------------------------------
C
      facts(0) = 1.0d0
      do i = 1,2*nterms
	 facts(i) = facts(i-1)*dsqrt(i+0.0d0)
      enddo
c
C     compute necessary powers of -i
C
      zeye(0) = 1.0D0
      do i = 1,nthmax
         zeye(i) = zeye(i-1)*dcmplx(0.0d0,-1.0d0)
      enddo
c
      do nm = 0,nterms
         do mth = -nm,nm
            localtmp(nm,mth) = 0.0d0
         enddo
      enddo
C
C     Loop over lambdas computing contributions to 
C     local expansion
C
      ntot = 1
      do nl = 1,nlambs
C
C     compute powers of lambda_nl
C
         rlampow(0) = whts(nl)
         rmul = rlams(nl)/(boxsize*scale)
         do nm = 1,nterms
            rlampow(nm) = rlampow(nm-1)*rmul
         enddo
c
C     add contributions to LOCAL expansion.
c
         do nm = 0,nterms
            mmax = numtets(nl)-1
            if (mmax.gt.nm) mmax = nm
	    rmul = rlampow(nm)
            do mth = 0,mmax
               ncurrent = ntot+mth
ccc               local(nm,mth) = local(nm,mth) + rmul*
               localtmp(nm,mth) = localtmp(nm,mth) + rmul*
     1          mexpdown(ncurrent)
            enddo
         enddo
         ntot = ntot+numtets(nl)
      enddo
c
C     scale the expansions according to formula
C
      do nm = 0,nterms
         do mth = 0,nm
            rscale = 1.0d0/(facts(nm-mth)*facts(nm+mth)) 
ccc            local(nm,mth) = local(nm,mth)*zeye(mth)*rscale/boxsize
            localtmp(nm,mth) = localtmp(nm,mth)*zeye(mth)*rscale/boxsize
         enddo
         do mth = 1,nm
            localtmp(nm,-mth) = dconjg(localtmp(nm,mth))
         enddo
      enddo
c
      do nm = 0,nterms
         do mth = -nm,nm
            local(nm,mth) = local(nm,mth) + localtmp(nm,mth)
         enddo
      enddo
      return
      end
