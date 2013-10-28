c***********************************************************************
      subroutine formmp_mindlinb(rscale,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,nparts,center,
     2           nterms,mpwork,chwork,quadvec,mpole)
c
c***********************************************************************
c     This subroutine forms a multipole expansion due to 
c     Mindlin B single and double layer sources.
c
c     INPUT:
c
c     rscale: scaling parameter for multipole expansion 
c             (typically the inverse boxsize)
c     rlame:  Lame coeffcients (rlam,rmu) 
c     source: source locations
c     ifsingle: SLP flag (1 means present)
c     sigma_sl: single layer force vector
c     ifdouble: DLP flag (1 means present)
c     dipstr:   DLP displacement jump
c     dipvec:   normal vector
c     nparts:   number of sources
c     center:   expansion center
c     nterms:   expansion length
c     mpwork:   multipole workspace 
c     chwork:   complex charge workspace 
c     quadvec:  quadrupole vector workspace
c
c     OUTPUT:
c
c     mpole:   multipole expansion 
c-----------------------------------------------------------------------
      implicit none
      integer i,j,ier
      real *8 rscale,rfac1,rr
      real *8 rlame(2)
      real *8 source(3,nparts)
      integer ifsingle
      real *8 sigma_sl(3,nparts)
      integer ifdouble
      real *8 dipstr(3,nparts)
      real *8 dipvec(3,nparts)
      integer nparts
      real *8 center(3)
      integer nterms
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpwork(0:nterms,-nterms:nterms)
      complex *16 chwork(nparts)
      real *8 quadvec(6,nparts)
c
      do i = 0,nterms
         do j = -nterms,nterms
            mpole(i,j) = 0.0d0
         enddo
      enddo
c
      if (ifsingle.eq.1) then
         do i = 1,nparts
            chwork(i) = 1.0d0/rscale
         enddo
         call l3dformmp_dp(ier,rscale,source,chwork,sigma_sl,nparts,
     1        center,nterms,mpwork)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole(i,j) = mpole(i,j) + mpwork(i,j)
            enddo
         enddo
      endif
c
      if (ifdouble.eq.1) then
         do i = 1,nparts
            rfac1 = 2.0d0*rlame(2)/rscale
            quadvec(1,i) = dipstr(1,i)*dipvec(1,i)
            quadvec(2,i) = dipstr(2,i)*dipvec(2,i)
            quadvec(3,i) = dipstr(3,i)*dipvec(3,i)
            quadvec(4,i) = (dipstr(1,i)*dipvec(2,i) +
     1                     dipstr(2,i)*dipvec(1,i))
            rr = quadvec(1,i)+quadvec(2,i)+quadvec(3,i)
            rr = rr/rscale
            quadvec(1,i) = quadvec(1,i)*rfac1
            quadvec(2,i) = quadvec(2,i)*rfac1
            quadvec(3,i) = -quadvec(3,i)*rfac1
            quadvec(3,i) = quadvec(3,i) - rlame(1)*rr*2.0d0
            quadvec(4,i) = quadvec(4,i)*rfac1
            quadvec(5,i) = 0.0d0
            quadvec(6,i) = 0.0d0
         enddo
         call l3dformmp_quad(ier,rscale,source,quadvec,nparts,
     1        center,nterms,mpwork)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole(i,j) = mpole(i,j) + mpwork(i,j)
            enddo
         enddo
      endif
      return
      end
c
c
c
c***********************************************************************
      subroutine formmp_mindlinb_trunc(rscale,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,nparts,center,
     2           nterms,mpwork,chwork,quadvec,mpole,wlege,nlege)
c
c***********************************************************************
c     This subroutine forms a multipole expansion due to 
c     Mindlin B single and double layer sources.
c
c     INPUT:
c
c     rscale: scaling parameter for multipole expansion 
c             (typically the inverse boxsize)
c     rlame:  Lame coeffcients (rlam,rmu) 
c     source: source locations
c     ifsingle: SLP flag (1 means present)
c     sigma_sl: single layer force vector
c     ifdouble: DLP flag (1 means present)
c     dipstr:   DLP displacement jump
c     dipvec:   normal vector
c     nparts:   number of sources
c     center:   expansion center
c     nterms:   expansion length
c     mpwork:   multipole workspace 
c     chwork:   complex charge workspace 
c     quadvec:  quadrupole vector workspace
c
c     OUTPUT:
c
c     mpole:   multipole expansion 
c-----------------------------------------------------------------------
      implicit none
      integer i,j,ier,nlege
      real *8 rscale,rfac1,rr
      real *8 rlame(2)
      real *8 source(3,nparts)
      integer ifsingle
      real *8 sigma_sl(3,nparts)
      integer ifdouble
      real *8 dipstr(3,nparts)
      real *8 dipvec(3,nparts)
      integer nparts
      real *8 center(3)
      integer nterms
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpwork(0:nterms,-nterms:nterms)
      complex *16 chwork(nparts)
      real *8 quadvec(6,nparts)
      real *8 wlege(*)
      real *8, allocatable :: dipvec3(:,:)
c
        if( nparts .le. 0 ) return
c
      do i = 0,nterms
         do j = -nterms,nterms
            mpole(i,j) = 0.0d0
         enddo
      enddo
c
      if (ifsingle.eq.1) then
         do i = 1,nparts
            chwork(i) = 1.0d0/rscale
         enddo
         call l3dformmp_dp_trunc
     $      (ier,rscale,source,chwork,sigma_sl,nparts,
     1       center,nterms,nterms,mpwork,wlege,nlege)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole(i,j) = mpole(i,j) + mpwork(i,j)
            enddo
         enddo
      endif
c
      if (ifdouble.eq.1) then
         do i = 1,nparts
            rfac1 = 2.0d0*rlame(2)/rscale
            quadvec(1,i) = dipstr(1,i)*dipvec(1,i)
            quadvec(2,i) = dipstr(2,i)*dipvec(2,i)
            quadvec(3,i) = dipstr(3,i)*dipvec(3,i)
            quadvec(4,i) = (dipstr(1,i)*dipvec(2,i) +
     1                     dipstr(2,i)*dipvec(1,i))
            rr = quadvec(1,i)+quadvec(2,i)+quadvec(3,i)
            rr = rr/rscale
            quadvec(1,i) = quadvec(1,i)*rfac1
            quadvec(2,i) = quadvec(2,i)*rfac1
            quadvec(3,i) = -quadvec(3,i)*rfac1
            quadvec(3,i) = quadvec(3,i) - rlame(1)*rr*2.0d0
            quadvec(4,i) = quadvec(4,i)*rfac1
            quadvec(5,i) = 0.0d0
            quadvec(6,i) = 0.0d0
         enddo
         call l3dformmp_quad_trunc(ier,rscale,source,quadvec,nparts,
     1        center,nterms,mpwork,wlege,nlege)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole(i,j) = mpole(i,j) + mpwork(i,j)
            enddo
         enddo
      endif
      return
      end
c
c
c
c***********************************************************************
      subroutine mtx_mindlinb(rscale,rlame,mpole,nterms,nlams,numfour,
     1           nexptot,mexpf,rlsc,nthmax,fexpe,nexpe,rlams,whts,
     2           numphys,nthmaxp,nexptotp,mexpphys)
c***********************************************************************
c     This subroutine converts a Mindlin B multipole expansion to
c     the corresponding plane wave expansion.
c
c     INPUT:
c
c     rscale:   multipole expansion scaling factor
c     rlame:    Lame coeffcients (rlam,rmu) 
c     mpole:    multipole expansion 
c     nterms:   expansion length
c     nlams:    number of pts in lambda integral quadrature
c     numfour:  number of azimuthal modes for each lambda (Fourier)
c     nexptot:  total number of Fourier modes
c     mexpf:    workspace to hold Fourier plane wave expansion
c     rlsc:     pecomputed array of coefficients needed by mpoletoexp
c     nthmax:   max number ofazimuthal (Fourier) modes
c     fexpe:    precomputed array of coefficients needed by ftophys
c     nexpe:    length of fexpe
c     rlams:    lmabda quadrature nodes
c     whts:     lmabda quadrature weights
c     numphys:  number of azimuthal modes for each lambda (physical)
c     nthmaxp:  max number ofazimuthal (physical) modes
c     nexptotp: total number of physical modes
c
c     OUTPUT:
c
c     mexpphys: outgoing plane wave expansion (physical)
c-----------------------------------------------------------------------
      implicit none
      real *8 rscale
      real *8 rlame(2)
      integer nterms
      complex *16 mpole(0:nterms,-nterms:nterms)
      integer nlams,numfour(nlams)
      integer nexptot
      complex *16 mexpf(nexptot)
      real *8  rlsc(nlams,0:nterms,0:nterms)
      integer nthmax,nthmaxp,nexptotp,nexpe
      complex *16 mexpphys(nexptotp)
      complex *16 fexpe(nexpe)
      real *8  rlams(nlams)
      real *8  whts(nlams),boxsize
      integer numphys(nlams)
c
      call mpoletoexp(mpole,nterms,nlams,numfour,nexptot,
     1        mexpf,rlsc)
      call ftophys(mexpf,nexptot,nlams,rlams,numfour,numphys,
     1             nthmax,mexpphys,nexptotp,fexpe,nexpe)
      boxsize = 1.0d0
      call mindlinbpwdiagall(rscale,rlame,boxsize,rlams,
     1        whts,nlams,numphys,nthmaxp,nexptotp,mexpphys)
c
      return
      end
c
c
c***********************************************************************
      subroutine mindlinbpwdiagall(rscale,rlame,boxsize,rlams,
     2           whts,nlams,numphys,nthmaxp,nexptotp,planewave)
c
c***********************************************************************
      implicit none
      real *8 rlame(2),rscale,pi,rlam,rmu,alpha,scfac,boxsize
      complex *16 planewave(nexptotp)
      real *8 rlams(nlams),whts(nlams)
      integer nlams,numphys(nlams)
      integer nthmaxp,nexptotp,ii,i,j
      complex *16 zfac,eye
c
c     After a multipole expansion for the full SLP force vector
c     has been converted to a plane wave expansion, that expansion
c     must be modified in diagonal form so that it then
c     represents the displacement. Unlike MINDLINBPWDIAG, all components
c     of displacement are available on output as gradient of plane wave
c     representation.
c     It doesn't include the 1/rmu scaling.
c
c     INPUT:
c
c     rlame(2)   Lame coefficients supplied in the form
c                rlame(1) = rlam, rlame(2) = rmu
c     rscale     scaling parameter
c     rlams      discretization points in lambda integral 
c     whts       discretization weights in lambda integral 
c     nlams      number of discretization pts. in lambda integral
c     numphys(j) number of Fourier modes needed in expansion
c                 expansion of alpha variable for lambda_j. 
c     nthmaxp    max_j numphys(j)
c     nexptotp   sum_j numphys(j)
c 
c     OUTPUT:
c                
c     planewave  planewave expansion
c-----------------------------------------------------------------------
c      
      pi = 4*datan(1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      zfac = 0.0d0
c
      rlam = rlame(1)
      rmu = rlame(2)
      alpha = (rlam+rmu)/(rlam+2*rmu)
      scfac = (1.0d0-alpha)/alpha/rscale
c
      ii = 1
      do i=1,nlams
      do j=1,numphys(i)
            zfac = 1.0d0/(rlams(i)**2)
            planewave(ii) = planewave(ii)*zfac*scfac
            ii = ii+1
      enddo
      enddo
      return
      end
c
