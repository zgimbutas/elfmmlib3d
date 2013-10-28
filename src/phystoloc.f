c
c     source to local translation operator using plane wave basis
c
c     The steps involved are:
c     
c     1) translate expansion in X,Y,Z directions,
c        calling subroutine XYZSHIFT
c
c     2) convert physical space representation to Fourier expansion,
c        calling subroutine PHYSTOF
c
c     3) convert Fourier representation to local expansion,
c        calling subroutine EXPTOLOCALD
c
c***********************************************************************
      subroutine mexpphystoloc(boxsize,x0y0z0,rlams,
     $        whts,nlams,numfour,numphys,nthmax,nexptot,
     $        nexptotp,mexpf,mexpphys,mexpnew,fexpback,nexpback,
     $        rlampow,facts,zeye,center,nterms,local,scloc)
c***********************************************************************
c
c     INPUT:
c
c     boxsize   box dimension in FMM structure
c     x0y0z0    old center
c     rlams     plane wave discretization nodes
c     whts      plane wave discretization weights
c     nlams     number of plane wave discretization points
c     numfour   number of Fourier modes in alpha integral  
c     numphys   number of <<physical>> samples of alpha integral  
c     nthmax    max of numfour
c     nexptot   tot number of modes in numfour
c     nexptotp  tot number of points in numphys
c     mexpf     workspace for Fourier exponential expansion
c     mexpphys  physical space exponential expansion
c     mexpnew   workspace for shifted physical space exp. expansion
c     fexpback  precomputed table for subroutine phystof
c     nexpback  length of fexpback
c     rlampow   work array
c     facts     real work array
c     zeye      complex work array
c     center    coordinates of center of local expansion
c
c     nterms    order of local expansion
c     scloc     scaling parameter for local expansion
c
c     OUTPUT:
c
c     local     local expansion about <<center>>
c
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nlams,nexptot,nexptotp,nthmax,nexpback
      integer numfour(nlams)
      integer numphys(nlams)
      real *8 x0y0z0(3)
      real *8 center(3)
      real *8 xdis,ydis,zdis,boxsize,scloc
      real *8 rlampow(0:nthmax)
      real *8 facts(0:2*nterms)
      complex *16 fexpback(nexpback)
      complex *16 mexpf(nexptot)
      complex *16 mexpnew(nexptotp)
      complex *16 mexpphys(nexptotp)
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 zeye(0:nthmax)
      real *8 rlams(nlams),whts(nlams)
c
      xdis = (center(1)-x0y0z0(1))/boxsize
      ydis = (center(2)-x0y0z0(2))/boxsize
      zdis = (center(3)-x0y0z0(3))/boxsize
      call xyzshift(mexpphys,rlams,nlams,numphys,mexpnew,nexptot,
     $             xdis,ydis,zdis)
      call phystof(mexpf,nexptot,nlams,rlams,numfour,numphys,nthmax,
     $             mexpnew,nexptotp,fexpback,nexpback)
      call exptolocald(local,nterms,rlams,whts,nlams,
     $             numfour,nthmax,nexptot,mexpf,boxsize,scloc,
     $             rlampow,facts,zeye)
      return
      end
