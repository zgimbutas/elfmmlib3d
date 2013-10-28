c***********************************************************************
      subroutine formmp_mindlinc(rscale,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,nparts,center,
     2           nterms,mpwork,chwork,quadvec,dipvec2,mpole1,mpole2)
c
c***********************************************************************
c     This subroutine forms a multipole expansion due to 
c     Mindlin C single and double layer sources.
c
c     INPUT:
c
c     rscale: scaling parameter for multipole expansion 
c             MUST BE INVERSE OF BOXSIZE IN FMM DATA STRUCTURE
c             IN ORDER FOR SUBSEQUENT CALL TO mtx_mindlinc.
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
c     dipvec:   dipole vector workspace
c
c     OUTPUT:
c
c     mpole1:   multipole expansion 
c     mpole2:   multipole expansion for correction to U_3
c-----------------------------------------------------------------------
      implicit none
      integer i,j,ier
      real *8 rscale,rfac1,rr,alpha
      real *8 rlame(2)
      real *8 source(3,nparts)
      integer ifsingle
      real *8 sigma_sl(3,nparts)
      integer ifdouble
      real *8 dipstr(3,nparts)
      real *8 dipvec(3,nparts)
      real *8 dipvec2(3,nparts)
      integer nparts
      real *8 center(3)
      integer nterms
      complex *16 mpole1(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
      complex *16 mpwork(0:nterms,-nterms:nterms)
      complex *16 chwork(nparts)
      real *8 quadvec(6,nparts)
c
      real *8, allocatable :: sigma_sl0(:,:)
c
      do i = 0,nterms
         do j = -nterms,nterms
            mpole1(i,j) = 0.0d0
            mpole2(i,j) = 0.0d0
         enddo
      enddo
c
      if( nparts .le. 0 ) return
c 
      alpha = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
c
      if (ifsingle.eq.1) then
         allocate( sigma_sl0(3,nparts) )
         do i = 1,nparts
            sigma_sl0(1,i) = +sigma_sl(1,i)
            sigma_sl0(2,i) = +sigma_sl(2,i)
            sigma_sl0(3,i) = -sigma_sl(3,i)
         enddo

         do i = 1,nparts
            chwork(i) = -alpha*source(3,i)
         enddo
         call l3dformmp_dp
     $      (ier,rscale,source,chwork,sigma_sl0,nparts,
     $      center,nterms,mpwork)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
         do i = 1,nparts
            chwork(i) = (2.0d0-alpha)*sigma_sl0(3,i)
         enddo
         call l3dformmp(ier,rscale,source,chwork,nparts,
     1        center,nterms,mpwork)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
         do i = 1,nparts
            chwork(i) = -(2.0d0-alpha)
         enddo
         call l3dformmp_dp
     $       (ier,rscale,source,chwork,sigma_sl0,nparts,
     1        center,nterms,mpole2)

      endif
c
      if (ifdouble.eq.1) then
         do i = 1,nparts
            rfac1 = 2.0d0*alpha*source(3,i)*rlame(2)
            quadvec(1,i) = -rfac1*dipstr(1,i)*dipvec(1,i)
            quadvec(2,i) = -rfac1*dipstr(2,i)*dipvec(2,i)
            quadvec(3,i) = -rfac1*dipstr(3,i)*dipvec(3,i)
            quadvec(4,i) = -rfac1*(dipstr(1,i)*dipvec(2,i) +
     1                     dipstr(2,i)*dipvec(1,i))
            quadvec(5,i) = rfac1*(dipstr(1,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(1,i))
            quadvec(6,i) = rfac1*(dipstr(2,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(2,i))
         enddo
         call l3dformmp_quad(ier,rscale,source,quadvec,nparts,
     1        center,nterms,mpwork)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
         do i = 1,nparts
            rfac1 = (2.0d0 - 2.0d0*alpha)*rlame(2)
            chwork(i) = 1.0d0
            dipvec2(1,i) = -rfac1*(dipstr(1,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(1,i))
            dipvec2(2,i) = -rfac1*(dipstr(2,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(2,i))
            dipvec2(3,i) = rfac1*(dipstr(3,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(3,i))
            dipvec2(3,i) = dipvec2(3,i) - rlame(1)*2.0d0*(alpha-1)*
     1                     (dipstr(1,i)*dipvec(1,i)+
     1                     dipstr(2,i)*dipvec(2,i)+
     1                     dipstr(3,i)*dipvec(3,i))
         enddo
         call l3dformmp_dp(ier,rscale,source,chwork,dipvec2,nparts,
     1        center,nterms,mpwork)
c
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
c     now z-component correction
c
         do i = 1,nparts
            rfac1 = -2.0d0*(2.0d0-alpha)*rlame(2)
            quadvec(1,i) = rfac1*dipstr(1,i)*dipvec(1,i)
            quadvec(2,i) = rfac1*dipstr(2,i)*dipvec(2,i)
            quadvec(4,i) = rfac1*(dipstr(1,i)*dipvec(2,i) +
     1                     dipstr(2,i)*dipvec(1,i))
            quadvec(5,i) = -rfac1*(dipstr(1,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(1,i))
            quadvec(6,i) = -rfac1*(dipstr(2,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(2,i))
            quadvec(3,i) = rfac1*dipstr(3,i)*dipvec(3,i)
         enddo
         call l3dformmp_quad(ier,rscale,source,quadvec,nparts,
     1        center,nterms,mpwork)
c
         do i = 0,nterms
            do j = -nterms,nterms
               mpole2(i,j) = mpole2(i,j) + mpwork(i,j)
            enddo
         enddo
c
      endif
      return
      end
c
c
c
c
c***********************************************************************
      subroutine formmp_mindlinc_trunc(rscale,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,nparts,center,
     2           nterms,mpwork,chwork,quadvec,dipvec2,mpole1,mpole2,
     $           wlege,nlege)
c
c***********************************************************************
c     This subroutine forms a multipole expansion due to 
c     Mindlin C single and double layer sources.
c
c     INPUT:
c
c     rscale: scaling parameter for multipole expansion 
c             MUST BE INVERSE OF BOXSIZE IN FMM DATA STRUCTURE
c             IN ORDER FOR SUBSEQUENT CALL TO mtx_mindlinc.
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
c     dipvec:   dipole vector workspace
c
c     OUTPUT:
c
c     mpole1:   multipole expansion 
c     mpole2:   multipole expansion for correction to U_3
c-----------------------------------------------------------------------
      implicit none
      integer i,j,ier
      real *8 rscale,rfac1,rr,alpha
      real *8 rlame(2)
      real *8 source(3,nparts)
      integer ifsingle
      real *8 sigma_sl(3,nparts)
      integer ifdouble
      real *8 dipstr(3,nparts)
      real *8 dipvec(3,nparts)
      real *8 dipvec2(3,nparts)
      integer nparts
      real *8 center(3)
      integer nterms
      complex *16 mpole1(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
      complex *16 mpwork(0:nterms,-nterms:nterms)
      complex *16 chwork(nparts)
      real *8 quadvec(6,nparts)
      integer nlege
      real *8 wlege(*)
      real *8, allocatable :: dipvec3(:,:)
      real *8, allocatable :: dipvec4(:,:)
c
      real *8, allocatable :: sigma_sl0(:,:)
c
      do i = 0,nterms
         do j = -nterms,nterms
            mpole1(i,j) = 0.0d0
            mpole2(i,j) = 0.0d0
         enddo
      enddo
c
      if( nparts .le. 0 ) return
c
      alpha = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
c
      if (ifsingle.eq.1) then

         allocate( sigma_sl0(3,nparts) )
         do i = 1,nparts
            sigma_sl0(1,i) = +sigma_sl(1,i)
            sigma_sl0(2,i) = +sigma_sl(2,i)
            sigma_sl0(3,i) = -sigma_sl(3,i)
         enddo

         do i = 1,nparts
            chwork(i) = -alpha*source(3,i)
         enddo
         call l3dformmp_dp_trunc
     $      (ier,rscale,source,chwork,sigma_sl0,nparts,
     $      center,nterms,nterms,mpwork,wlege,nlege)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
         do i = 1,nparts
            chwork(i) = (2.0d0-alpha)*sigma_sl0(3,i)
         enddo
         call l3dformmp_trunc(ier,rscale,source,chwork,nparts,
     1        center,nterms,nterms,mpwork,wlege,nlege)
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
         do i = 1,nparts
            chwork(i) = -(2.0d0-alpha)
         enddo
         call l3dformmp_dp_trunc
     $       (ier,rscale,source,chwork,sigma_sl0,nparts,
     1        center,nterms,nterms,mpole2,wlege,nlege)

      endif
c
      if (ifdouble.eq.1) then
         do i = 1,nparts
            rfac1 = 2.0d0*alpha*source(3,i)*rlame(2)
            quadvec(1,i) = -rfac1*dipstr(1,i)*dipvec(1,i)
            quadvec(2,i) = -rfac1*dipstr(2,i)*dipvec(2,i)
            quadvec(3,i) = -rfac1*dipstr(3,i)*dipvec(3,i)
            quadvec(4,i) = -rfac1*(dipstr(1,i)*dipvec(2,i) +
     1                     dipstr(2,i)*dipvec(1,i))
            quadvec(5,i) = rfac1*(dipstr(1,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(1,i))
            quadvec(6,i) = rfac1*(dipstr(2,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(2,i))
         enddo
         call l3dformmp_quad_trunc(ier,rscale,source,quadvec,nparts,
     1        center,nterms,mpwork,wlege,nlege)

         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
         do i = 1,nparts
            rfac1 = (2.0d0 - 2.0d0*alpha)*rlame(2)
            chwork(i) = 1.0d0
            dipvec2(1,i) = -rfac1*(dipstr(1,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(1,i))
            dipvec2(2,i) = -rfac1*(dipstr(2,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(2,i))
            dipvec2(3,i) = rfac1*(dipstr(3,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(3,i))
            dipvec2(3,i) = dipvec2(3,i) - rlame(1)*2.0d0*(alpha-1)*
     1                     (dipstr(1,i)*dipvec(1,i)+
     1                     dipstr(2,i)*dipvec(2,i)+
     1                     dipstr(3,i)*dipvec(3,i))
         enddo
c         call l3dformmp_dp_trunc
c     $      (ier,rscale,source,chwork,dipvec2,nparts,
c     1        center,nterms,nterms,mpwork,wlege,nlege)
c       
         call l3dformmp_dipole_trunc
     $      (ier,rscale,source,dipvec2,nparts,
     1        center,nterms,mpwork,wlege,nlege)
c
         do i = 0,nterms
            do j = -nterms,nterms
               mpole1(i,j) = mpole1(i,j) + mpwork(i,j)
            enddo
         enddo
c
c     now z-component correction
c
         do i = 1,nparts
            rfac1 = 2.0d0*(2.0d0-alpha)*rlame(2)
            quadvec(1,i) = -rfac1*dipstr(1,i)*dipvec(1,i)
            quadvec(2,i) = -rfac1*dipstr(2,i)*dipvec(2,i)
            quadvec(3,i) = -rfac1*dipstr(3,i)*dipvec(3,i)
            quadvec(4,i) = -rfac1*(dipstr(1,i)*dipvec(2,i) +
     1                     dipstr(2,i)*dipvec(1,i))
            quadvec(5,i) = rfac1*(dipstr(1,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(1,i))
            quadvec(6,i) = rfac1*(dipstr(2,i)*dipvec(3,i) +
     1                     dipstr(3,i)*dipvec(2,i))
         enddo
         call l3dformmp_quad_trunc(ier,rscale,source,quadvec,nparts,
     1        center,nterms,mpwork,wlege,nlege)
c       
         do i = 0,nterms
            do j = -nterms,nterms
               mpole2(i,j) = mpole2(i,j) + mpwork(i,j)
            enddo
         enddo
c
      endif
      return
      end
c
c
c
c
c***********************************************************************
      subroutine mtx_mindlinc(rlame,mpole1,mpole2,nterms,nlams,numfour,
     1           nexptot,mexpf,rlsc,nthmax,fexpe,nexpe,rlams,whts,
     2           numphys,nthmaxp,nexptotp,mexpphys1,mexpphys2)
c***********************************************************************
c     This subroutine converts a Mindlin B multipole expansion to
c     the corresponding plane wave expansion.
c
c     INPUT:
c
c     rlame:    Lame coeffcients (rlam,rmu) 
c     mpole1:   main  multipole expansion 
c     mpole2:   correction for U_3 multipole expansion 
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
c     mexpphys1: main outgoing plane wave expansion (physical)
c     mexpphys2: correction outgoing plane wave expansion (physical)
c-----------------------------------------------------------------------
      implicit none
      real *8 rlame(2)
      integer nterms
      complex *16 mpole1(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
      integer nlams,numfour(nlams)
      integer nexptot
      complex *16 mexpf(nexptot)
      real *8  rlsc(nlams,0:nterms,0:nterms)
      integer nthmax,nthmaxp,nexptotp,nexpe
      complex *16 mexpphys1(nexptotp)
      complex *16 mexpphys2(nexptotp)
      complex *16 fexpe(nexpe)
      real *8  rlams(nlams)
      real *8  whts(nlams)
      integer numphys(nlams)
c
      call mpoletoexp(mpole1,nterms,nlams,numfour,nexptot,
     1        mexpf,rlsc)
      call ftophys(mexpf,nexptot,nlams,rlams,numfour,numphys,
     1             nthmax,mexpphys1,nexptotp,fexpe,nexpe)
c
      call mpoletoexp(mpole2,nterms,nlams,numfour,nexptot,
     1        mexpf,rlsc)
      call ftophys(mexpf,nexptot,nlams,rlams,numfour,numphys,
     1             nthmax,mexpphys2,nexptotp,fexpe,nexpe)
c
      return
      end
c




