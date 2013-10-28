c**********************************************************************
      subroutine mindlinmpevalb(rlame,rscale,center,mpole,nterms,
     1		ztarg,skels,chargeskels,nskel,fld,
     2          ifhess,gradfld)
c**********************************************************************
c     This subroutine evaluates the displacement and strain
c     due to a Mindlin B multipole expansion.
c
c     fld = -gradient
c     gradfld = Hessian with gradfld(i,j) = dx_i( fld(j)).
c
c     where rscale defines scaling parameter.     
c
c     Subroutine for computing the displacement and gradient of 
c     displacement for a Mindlin B multipole expansion.
c     Since the multipole expansion itself requires the 
c     incorporation of two z-antiderivatives, the plane wave 
c     representation is needed as an intermediary step.
c     Thus, the algorithm is: 
c 
c     mpole -> plane wave  -> incorporate anti-derivatives -> 
c     evaluate all derivatives from plane-wave representation.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rlame  =    Lame coefficients
c     rscale =    scaling parameter 
c     center =    expansion center
c     mpole  =    multipole expansion in 2d matrix format
c     nterms =    order of the multipole expansion
c     ztarg  =    target location
C     skels     =  proxy charge locations
C     chargeskels  =  proxy charge strengths
C     nskel     =  number of proxy charges
C     ifhess    =  flag to compute gradient of displacement
C
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     fld       =    -gradient at ztarg
c     gradfld   =    gradient of fld (if requested)
c-----------------------------------------------------------------------
      implicit none
      integer ii,i,j,nterms,nlams,nthmax,nthmaxp,nskel
      integer ifhess,ier,iffld
      real *8 center(3),ztarg(3)
      real *8 rss
      real *8 rlame(2)
      real *8 skels(3,nskel)
      real *8 rscale,dx,dy,dz,rho,boxsize,rscale2,zfac,alpha,scfac
      real *8 boxsizem
      complex *16 chargeskels(nskel)
      complex *16 fld(3),cfield(3)
      complex *16 gradfld(3,3),chess(6),cpot
      complex *16 mpole(0:nterms,-nterms:nterms)
c
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      scfac  = (1.0d0-alpha)/alpha
      call  evalintkerbsc(skels,chargeskels,nskel,rlame,
     1     rscale,ztarg,fld,ifhess,gradfld)
      iffld = 1
      call l3dmpevalhess(rscale,center,mpole,nterms,ztarg,cpot,
     1         iffld,cfield,ifhess,chess,ier)
      rss = scfac/rscale
      fld(1) = fld(1) -rss*cfield(1)
      fld(2) = fld(2) -rss*cfield(2)
      fld(3) = fld(3) +rss*cfield(3)
      if (ifhess.eq.1) then
         gradfld(1,1) = gradfld(1,1) -rss*chess(1)
         gradfld(2,1) = gradfld(2,1) -rss*chess(4)
         gradfld(3,1) = gradfld(3,1) -rss*chess(5)
         gradfld(1,2) = gradfld(1,2) -rss*chess(4)
         gradfld(2,2) = gradfld(2,2) -rss*chess(2)
         gradfld(3,2) = gradfld(3,2) -rss*chess(6)
         gradfld(1,3) = gradfld(1,3) +rss*chess(5)
         gradfld(2,3) = gradfld(2,3) +rss*chess(6)
         gradfld(3,3) = gradfld(3,3) +rss*chess(3)
      endif
      return
      end
c
C**********************************************************************
      subroutine mpdiff(mpole,mpole2,nterms)
C
C     set  mpole = mpole2 - mpole
C-----------------------------------------------------------------------
C     INPUT:
C
C     mpole   =    multipole expansion 
C     mpole2  =    multipole expansion to be subtracted
C     nterms  =    order of the multipole expansion
C
c     OUTPUT:
C
C     mpole  =    mpole2-mpole
C
C-----------------------------------------------------------------------
      implicit none
      integer i,j,nterms
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
c
      do i = 0,nterms
         do j = -i,i
            mpole(i,j) = mpole2(i,j) - mpole(i,j)
         enddo
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine mindlinmpevalc(rlame,rscale,center,mpole,mpole2,
     1		nterms,ztarg,fld,ifhess,gradfld)
c**********************************************************************
c     This subroutine evaluates the gradient of the 
c     multipole expansion due to Mindlin C sources.
c
c     fld = -gradient
c     gradfld = Hessian with gradfld(i,j) = dx_i( fld(j)).
c
c     where rscale defines scaling parameter.     
c
c     Subroutine for computing the displacement and gradient of 
c     displacement for Mindlin C multipole expansions.
c-----------------------------------------------------------------------
c     INPUT:
c
c     rlame  =    Lame coefficients
c     rscale =    scaling parameter 
c     center =    expansion center
c     mpole  =    first C multipole expansion in 2d matrix format
c     mpole2 =    second C multipole expansion in 2d matrix format
c     nterms =    order of the multipole expansion
c     ztarg  =    target location
C     ifhess    =  flag to compute gradient of displacement
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     fld       =    -gradient at ztarg
c     gradfld   =    gradient of fld (if requested)
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ifhess,ier,iffld
      real *8 center(3),ztarg(3)
      real *8 rlame(2)
      real *8 rscale
      complex *16 fld(3)
      complex *16 cfield(3),chess(6),cfield2(3),cpot,cpot2
      complex *16 gradfld(3,3)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
c
      iffld = 1
      call l3dmpevalhess(rscale,center,mpole,nterms,ztarg,cpot,
     1         iffld,cfield,ifhess,chess,ier)
      call l3dmpeval(rscale,center,mpole2,nterms,ztarg,cpot2,
     1               iffld,cfield2,ier)
      fld(1) = -cfield(1)
      fld(2) = -cfield(2)
      fld(3) = -cfield(3)
      fld(3) = fld(3) -cpot2
      if (ifhess.eq.1) then
         gradfld(1,1) = -chess(1)
         gradfld(2,1) = -chess(4)
         gradfld(3,1) = -chess(5)
         gradfld(1,2) = -chess(4)
         gradfld(2,2) = -chess(2)
         gradfld(3,2) = -chess(6)
         gradfld(1,3) = -chess(5)
         gradfld(2,3) = -chess(6)
         gradfld(3,3) = -chess(3)
c
         gradfld(1,3) = gradfld(1,3) - cfield2(1) 
         gradfld(2,3) = gradfld(2,3) - cfield2(2) 
         gradfld(3,3) = gradfld(3,3) - cfield2(3) 
      endif
      return
      end
C
C
C
C
C
C***********************************************************************
      subroutine evalintkerbsc(source,charge1,ns,rlame,
     1     rscale,targ,fld,ifhess,gradfld)
C***********************************************************************
C
C     This subroutine computes a harmonic function defined by 
C     a set of undifferentiated scalar Mindlin B sources
C     generated by subroutine mptoslp.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     source    = source locations
C     charge1    = single layer force vector
C     ns        = number of sources
C     rlame     = Lame coefficients
C     rscale    = scaling parameter of original mpole expansion
C     targ      = target
C     ifhess    = flag to compute gradient of displacement
C                 1 means compute
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     fld       = displacement
C     gradlfd   = gradient of displacement (if requested)
C
C---------------------------------------------------------------------
      implicit none
      integer ns,i,ifhess
      real *8 source(3,ns),targ(3)
      real *8 rlame(2)
      real *8 rscale,r1,r2,r3,rr,rkerb,scfac,alpha,dd
      real *8 u1,u2,u3,u11,u12,u13,u22,u23,u33
      complex *16 charge1(ns)
      complex *16 fld(3), gradfld(3,3)
C
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      scfac  = (1.0d0-alpha)/alpha
      fld(1) = 0.0d0
      fld(2) = 0.0d0
      fld(3) = 0.0d0
      if (ifhess.eq.1) then
         gradfld(1,1) = 0.0d0
         gradfld(1,2) = 0.0d0
         gradfld(1,3) = 0.0d0
         gradfld(2,2) = 0.0d0
         gradfld(2,3) = 0.0d0
         gradfld(3,3) = 0.0d0
      endif
      do i = 1,ns
         r1 = targ(1)-source(1,i)
         r2 = targ(2)-source(2,i)
         r3 = targ(3)-source(3,i)
         rr = sqrt(r1*r1+r2*r2+r3*r3)
         dd = 1.0d0/(rr-r3)
         u1 = -r1*dd
         u2 = -r2*dd
         u3 = dlog(rr-r3)
         u11 = dd*(-1.0d0+r1*r1*dd/rr)
         u12 = r1*r2*dd*dd/rr
         u13 = r1*dd/rr
         u22 = dd*(-1.0d0+r2*r2*dd/rr)
         u23 = r2*dd/rr
         u33 = 1.0d0/rr
         fld(1) = fld(1) + u1*charge1(i)
         fld(2) = fld(2) + u2*charge1(i)
         fld(3) = fld(3) + u3*charge1(i)
         if (ifhess.eq.1) then
            gradfld(1,1) = gradfld(1,1) - u11*charge1(i)
            gradfld(1,2) = gradfld(1,2) - u12*charge1(i)
            gradfld(1,3) = gradfld(1,3) - u13*charge1(i)
            gradfld(2,2) = gradfld(2,2) - u22*charge1(i)
            gradfld(2,3) = gradfld(2,3) - u23*charge1(i)
            gradfld(3,3) = gradfld(3,3) - u33*charge1(i)
         endif
      enddo
      fld(1) = fld(1)*scfac*rscale
      fld(2) = fld(2)*scfac*rscale
      fld(3) = fld(3)*scfac*rscale
      if (ifhess.eq.1) then
         gradfld(1,1) = gradfld(1,1)*scfac*rscale
         gradfld(1,2) = gradfld(1,2)*scfac*rscale
         gradfld(1,3) = gradfld(1,3)*scfac*rscale
         gradfld(2,1) = gradfld(1,2)
         gradfld(2,2) = gradfld(2,2)*scfac*rscale
         gradfld(2,3) = gradfld(2,3)*scfac*rscale
         gradfld(3,1) = -gradfld(1,3)
         gradfld(3,2) = -gradfld(2,3)
         gradfld(3,3) = -gradfld(3,3)*scfac*rscale
      endif
      return
      end
C
c
C**********************************************************************
      subroutine mpdzminus2(mpole,nterms)
C
C     convert mpole to dz^{-2} mpole, assuming appropriate border
C     of mpole is identically zero.
C
C-----------------------------------------------------------------------
C     INPUT:
C
C     mpole   =    multipole expansion 
C     nterms =    order of the multipole expansion
C
c     OUTPUT:
C
C     mpole  =    dz^{-2} mpole
C
C-----------------------------------------------------------------------
      implicit none
      integer i,j,nterms
      real *8 rs
      complex *16 mpole(0:nterms,-nterms:nterms)
c
      do j = -nterms+2,nterms-2
         do i = abs(j)+2,nterms
            rs = 1.0d0/sqrt((i-j)*(i-j-1.0d0)*(i+j)*(i+j-1.0d0))
            mpole(i-2,j) = mpole(i,j)*rs
         enddo
      enddo
      do i = nterms-1,nterms
      do j = -i,i
         mpole(i,j) = 0.0d0
      enddo
      enddo
      return
      end
c
c**********************************************************************
      subroutine mindlinbmploc(rscale,center,mpole,nterms,radius,
     1		rscale2,center2,local,nterms2,skels,chargeskels,nskel,
     2          lwork,rlame)
c**********************************************************************
c     This subroutine converts a Mindlin B outgoing representation
c     (mpole+skels/chargeskels to a local expansion
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale      =    scaling parameter for mpole
c     center      =    expansion center
c     mpole       =    multipole expansion component of Mindlin B
c     nterms      =    order of the multipole expansion
c     radius      =    radius for use in constructing local expansion
c                      from skel sources
c     rscale2     =    scaling parameter for local expansion
c     center2     =    expansion center
C     skels       =    proxy charge locations
C     chargeskels =    proxy charge strengths
C     nskel       =    number of proxy charges
c     lwork       =    workspace for a local expansion
c     rlame       =    Lame coefficients
c
c     OUTPUT:
c
c     local       =    local expansion for Mindlin B
c     nterms2     =    order of the local expansion
c
c
      implicit none
      integer nterms,nterms2,nskel,ier,i,j
      real *8 center(3),center2(3),rlame(2)
      real *8 skels(3,nskel)
      real *8 radius,rscale,rscale2,alpha,scfac
      complex *16 chargeskels(nskel)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 lwork(0:nterms2,-nterms2:nterms2)
c
      call l3dformtaintkerbf90sc(ier,rscale2,rscale,skels,
     $            chargeskels,nskel,rlame,center2,
     $            radius,nterms2,local)
c
      call l3dmplocquadu(rscale,center,mpole,nterms,
     1            rscale2,center2,lwork,nterms2,ier)
    
      alpha  = (rlame(1)+rlame(2))/(rlame(1)+2*rlame(2))
      scfac  = (1.0d0-alpha)/alpha/rscale
c
      do i = 0,nterms2
      do j = -i,i   
         local(i,j) = local(i,j) + lwork(i,j)*scfac
      enddo
      enddo
      return
      end
