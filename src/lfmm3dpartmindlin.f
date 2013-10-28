cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating the Mindlin B,C image 
c     contributions.  (FORTRAN 90 VERSION)
c
c     lfmm3dmindlinparttarg - Laplace FMM in R^3: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     l3dmindlinpartdirect - direct interactions in R^3: evaluate interactions
c       source to target via direct O(N^2) algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the Laplace particle FMM in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
      subroutine lfmm3dmindlinparttarg(ier,iprec,nsource,source,
     $     rlame,ifsingle,sigma_sl,ifdouble,dipstr,dipvec,
     $     ntarget,target,ifptfrctarg,ptfrctarg,ifhesstarg,hesstarg)
c       
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dparttargmain.
c
c       INPUT:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c
c       nsource     :   number of sources
c       source      :   source locations
c       rlame       :   Lame coefficients
c       ifsingle    :   single layer (SLP) FLAG
c                       ifsingle = 1   =>  include SLP
c                                     otherwise do not
c       sigma_sl    :   SLP vector
c       ifdouble    :    double layer (DLP) FLAG
c                       ifdouble = 1   =>  include DLP 
c                                     otherwise do not
c       dipstr      :   DLP strength vector
c       dipvec      :   dipole orientation vectors
c       ntarget     :   number of targets
c       target      :   target locations
c       ifptfrctarg :   target displacement flag 
c                       (1=compute ptfrctarg, otherwise no)
c       ifhesstarg  :   target strain flag 
c                       (1=compute strain, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       ptfrctarg   : displacement at target locations 
c       hesstarg    : strain at target locations 
c
      implicit none
      integer ier,lused7,iprec,nsource,ifsingle,ifdouble,ntarget
      integer ifptfrctarg,ifhesstarg
      real *8 rlame(2)
      real *8 source(3,nsource)
      real *8 sigma_sl(3,nsource)
      real *8 dipstr(3,nsource)
      real *8 dipvec(3,nsource)
      real *8 target(3,nsource)
      complex *16 ptfrctarg(3,ntarget)
      complex *16 hesstarg(3,3,ntarget)
c
      real *8 timeinfo(10)
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 evels of refinment, which is 
c     more than enough for any non-pathological case.
c
 
      integer laddr(2,200)
      integer nterms(0:200)
      integer box(20)
      integer box1(20)
      integer ichargesort,idipstrsort,idipvecsort,i
      integer ifevalfar,ifevalloc,ifprint,ihesstarg,ii,iisource
      integer iitarget,imptemp,imptemp2,iptfrctarg,irmlexp
      integer isourcesort,itargetsort,iwlists,lchargesort,ldipstr
      integer iprec1,ldipvec,lhesstarg,lmptemp,lmptot,lptfrctarg
      integer lsourcesort,ltargetsort,lused,lwlists,nbox
      integer nboxes,nchwork,nexpback,nexpe,nexptot,nexptotp
      integer nlams,nlev,nmax,nskel,nthmax,nthmaxp,ntot
      real *8 boxsize,done,epsfmm,pi,rsize
      real *8 bsize(0:200)
      real *8 scale(0:200)
      real *8 center(3)
      real *8 center0(3),corners0(3,8)
      real *8 center1(3),corners1(3,8)
      complex *16 ima
      complex *16 ptemp,ftmp(3)
c
      integer, allocatable :: numfour(:)
      integer, allocatable :: numphys(:)
      integer, allocatable :: iaddr(:)
      real *8, allocatable :: rlams(:)
      real *8, allocatable :: whts(:)
      real *8, allocatable :: w(:)
      real *8, allocatable :: wlists(:)
      real *8, allocatable :: wrmlexp(:)
      real *8, allocatable :: rlsc(:,:,:)
      real *8, allocatable :: vec3(:,:)
      real *8, allocatable :: quadvec(:,:)
      real *8, allocatable :: skels(:,:,:)
      complex *16, allocatable :: chwork(:)
      complex *16, allocatable :: fexpe(:)
      complex *16, allocatable :: fexpback(:)
      complex *16, allocatable :: mexpf(:)
      complex *16, allocatable :: mexpphys(:)
      complex *16, allocatable :: chargeskel(:,:)
      complex *16, allocatable :: cw(:)
c       
      data ima/(0.0d0,1.0d0)/
c       
      ier=0
      lused7 = 0
c       
      done=1
      pi=4*atan(done)
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
      ifprint=1
c
c     set fmm tolerance based on iprec flag.
c
      if( iprec .eq. -2 ) epsfmm=.5d-0 
      if( iprec .eq. -1 ) epsfmm=.5d-1
      if( iprec .eq. 0 ) epsfmm=.5d-2
      if( iprec .eq. 1 ) epsfmm=.5d-3
      if( iprec .eq. 2 ) epsfmm=.5d-6
      if( iprec .eq. 3 ) epsfmm=.5d-9
      if( iprec .eq. 4 ) epsfmm=.5d-12
      if( iprec .eq. 5 ) epsfmm=.5d-15
      if( iprec .eq. 6 ) epsfmm=0
c       
      if (ifprint .ge. 1) call prin2('epsfmm=*',epsfmm,1)
c
c
c     set criterion for box subdivision (number of sources per box)
c
      if( iprec .eq. -2 ) nbox=40
      if( iprec .eq. -1 ) nbox=50
      if( iprec .eq. 0 ) nbox=80
      if( iprec .eq. 1 ) nbox=160
      if( iprec .eq. 2 ) nbox=400
      if( iprec .eq. 3 ) nbox=800
      if( iprec .eq. 4 ) nbox=1200
      if( iprec .eq. 5 ) nbox=1400
      if( iprec .eq. 6 ) nbox=nsource+ntarget
      nbox = 40
c
      if (ifprint .ge. 1) call prinf('nbox=*',nbox,1)
c
c
c     create oct-tree data structure
c
      ntot = 100*(nsource+ntarget)+10000
      do ii = 1,10
         allocate (wlists(ntot))
         call lfmm3dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,rsize,
     $        wlists,ntot,lused7)
         if (ier.ne.0) then
            deallocate(wlists)
            ntot = ntot*1.5
            call prinf(' increasing allocation, ntot is *',ntot,1)
         else
           exit
         endif
      enddo
ccc      call prin2(' center is *',center,3)
ccc      call prin2(' rsize is *',rsize,1)
      if (ier.ne.0) then
         call prinf(' exceeded max allocation, ntot is *',ntot,1)
         ier = 4
         return
      endif
c
      do i = 0,nlev
         boxsize=rsize/2.0d0**i
         scale(i) = 1.0d0/boxsize
      enddo
      if (ifprint .ge. 1) call prin2('scale=*',scale,nlev+1)
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
c     carve up real workspace 
c
c     isourcesort is pointer for sorted source coordinates
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
      lused7=1
      isourcesort = 1
      lsourcesort = 3*nsource
      itargetsort = isourcesort+lsourcesort
      ltargetsort = 3*ntarget
      ichargesort = itargetsort+ltargetsort
      if (ifsingle.eq.1) then
         lchargesort = 3*nsource
      else
         lchargesort = 3
      endif
      idipvecsort = ichargesort+lchargesort
      if (ifdouble.eq.1) then
         ldipvec = 3*nsource
         ldipstr = 3*nsource
      else
         ldipvec = 3
         ldipstr = 3
      endif
      idipstrsort = idipvecsort + ldipvec
      lused7 = idipstrsort + ldipstr
      allocate(w(lused7),stat=ier)
      if (ier.ne.0) then
         call prinf(' cannot allocate bulk real FMM workspace,
     1                  lused7 is *',lused7,1)
         ier = 8
         return
      endif
c
c       ... allocate the potential and field arrays
c
      iptfrctarg = 1
      lptfrctarg = 3*ntarget
      ihesstarg=iptfrctarg+lptfrctarg
      if( ifhesstarg .eq. 1) then
        lhesstarg = (3*3*ntarget)
      else
        lhesstarg=9
      endif
      lused7=ihesstarg+lhesstarg
c      
      if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c      
      nmax = 0
      do i = 0,nlev
         bsize(i)=rsize/2.0d0**i
         call l3dterms(epsfmm, nterms(i), ier)
         if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
      enddo
c
c     work arrays for forming expansions
c
      nchwork = max(nbox,4*nmax*nmax)
      allocate(vec3(3,nchwork))
      allocate(quadvec(6,nchwork))
      allocate(chwork(nchwork))
c
c     work arrays for plane wave expansions
c
      allocate(rlams(100))
      allocate(whts(100))
      iprec1 = iprec
      if (iprec.lt.0) iprec1=0
      call lwtsexp3(iprec1,rlams,whts,nlams)
      allocate(numfour(nlams))
      allocate(numphys(nlams))
      call numthetahalf(numfour,nlams)
      call numthetafour(numphys,nlams)
c
      do i=1,nlams
         numfour(i) = numfour(i)+2
         numphys(i) = numphys(i)+2
      enddo
      nexptot = 0
      nexptotp = 0
      nthmax = 0
      nthmaxp = 0
      do i=1,nlams
         nexptot = nexptot + numfour(i)
         if (numfour(i).gt.nthmax) nthmax = numfour(i)
           nexptotp = nexptotp + numphys(i)
         if (numphys(i).gt.nthmaxp) nthmaxp = numphys(i)
      enddo
      call getfexplengths(nlams,numfour,numphys,nexpe,nexpback)
      allocate(fexpe(nexpe))
      allocate(fexpback(nexpback))
      allocate(mexpf(nexptot))
      allocate(mexpphys(nexptotp))
      call mkfexp(nlams,numfour,numphys,fexpe,nexpe,
     1       fexpback,nexpback)
      allocate(rlsc(nlams,0:nmax,0:nmax))
      call rlscini(rlsc,nlams,rlams,nmax)
c
      if (ifprint .ge. 1) call prinf('nterms=*',nterms,nlev+1)
      if (ifprint .ge. 1) call prinf('nlams=*',nlams,1)
      if (ifprint .ge. 1) call prinf('nmax=*',nmax,1)
c
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(6,nboxes).
c
c     imptemp is pointer for single expansion (dimensioned by nmax)
c   
c       ... allocate temporary arrays
c
      imptemp = lused7
      lmptemp = (nmax+1)*(2*nmax+1)
      imptemp2 = imptemp+lmptemp
      lused7 = imptemp2+lmptemp
      allocate(cw(lused7),stat=ier)
      allocate(iaddr(6*nboxes),stat=ier)
      if (ier.ne.0) then
         call prinf(' cannot allocate bulk FMM workspace,
     1                  lused7 is *',lused7,1)
         ier = 8
         return
      endif
      nskel = 4*nmax
      allocate(skels(3,nskel,nboxes),stat=ier)
      allocate(chargeskel(nskel,nboxes),stat=ier)
c
c     reorder sources, targets so that each box holds
c     contiguous list of source/target numbers.
c
      call l3dpartreorderv2(nsource,source,ifsingle,sigma_sl,
     $     wlists(iisource),ifdouble,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 
c       
      call l3dreordertarg(ntarget,target,wlists(iitarget),
     1       w(itargetsort))
c
      if (ifprint .ge. 1) call prinf('finished reordering=*',ier,1)
      if (ifprint .ge. 1) call prinf('ier=*',ier,1)
      if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
      if (ifprint .ge. 1) call prinf('nlev=*',nlev,1)
      if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
      if (ifprint .ge. 1) call prinf('lused7=*',lused7,1)
c

c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call l3dpartmindlinalloc(wlists(iwlists),iaddr,nboxes,lmptot,
     1       nterms,nexptotp)
c
      if (ifprint .ge. 1) call prinf(' lmptot is *',lmptot,1)
c       
      irmlexp = 1
      lused7 = irmlexp + lmptot 
      if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
      allocate(wrmlexp(lused7),stat=ier)
      if (ier.ne.0) then
         call prinf(' cannot allocate mpole expansion workspace,
     1                  lused7 is *',lused7,1)
         ier = 16
         return
      endif
c
c       
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far field and local fields (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
      ifevalfar=1
      ifevalloc=1
c
      call lfmm3dmindlinparttargmain(ier,iprec,
     $   ifevalfar,ifevalloc,nsource,w(isourcesort),wlists(iisource),
     $   rlame,ifsingle,w(ichargesort),ifdouble,w(idipstrsort),
     $   w(idipvecsort),ntarget,
     $   w(itargetsort),wlists(iitarget),ifptfrctarg,cw(iptfrctarg),
     $   ifhesstarg,cw(ihesstarg),epsfmm,iaddr,wrmlexp(irmlexp),
     $   vec3,chwork,quadvec,nchwork,skels,chargeskel,nskel,
     $   cw(imptemp),cw(imptemp2),lmptemp,nboxes,laddr,nlev,
     $   scale,bsize,nterms,wlists(iwlists),lwlists,mexpf,nexptot,
     $   mexpphys,nexptotp,nlams,rlams,whts,numfour,numphys,fexpe,
     $   nexpe,fexpback,nexpback,rlsc,nmax)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
      if( ier .ne. 0 ) ier = 0
c
      if (ifprint .ge. 1) call prinf('lwlists=*',lused,1)
      if (ifprint .ge. 1) call prinf('lused total =*',lused7,1)
c       
      if (ifprint .ge. 1) 
     $    call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
      if(ifptfrctarg .eq. 1 )
     $   call l3dfsort(ntarget,wlists(iitarget),cw(iptfrctarg),
     $          ptfrctarg)
      if(ifhesstarg .eq. 1) 
     $   call l3dssort(ntarget,wlists(iitarget),cw(ihesstarg),
     $          hesstarg)
c       
      return
      end
c
        subroutine lfmm3dmindlinparttargmain(ier,iprec,
     $     ifevalfar,ifevalloc,nsource,sourcesort,isource,
     $     rlame,ifsingle,chargesort,ifdouble,dipstrsort,
     $     dipvecsort,ntarget,
     $     targetsort,itarget,ifptfrctarg,ptfrctarg,
     $     ifhesstarg,hesstarg,epsfmm,iaddr,rmlexp,
     $     vec3,chwork,quadvec,nchwork,skels,chargeskel,nskel,
     $     mptemp,mptemp2,lmptemp,nboxes,laddr,nlev,
     $     scale,bsize,nterms,wlists,lwlists,mexpf,nexptot,
     $     mexpphys,nexptotp,nlams,rlams,whts,numfour,numphys,fexpe,
     $     nexpe,fexpback,nexpback,rlsc,nmax)
C
        implicit none
        integer ier,iprec,ifevalfar,ifevalloc,nsource
        integer ifsingle,ifdouble,ntarget,ifptfrctarg,lwlists
        integer ifhesstarg,nchwork,nskel,lmptemp,nboxes,nlev
        integer nexptot,nexptotp,nlams,nexpe,nexpback
        integer numfour(nlams),numphys(nlams)
        integer nmax,i,ibox,ifdirect2,ifdirect3,ifdirect4
        integer ifprint,ilev,ilist,itype,j,jbox,jj,k,kk,level1,ll
        integer nkids,npts,level,nlist,nops,nquad,nthmax,nthmaxp
        integer iaddr(6,nboxes)
        integer isource(nsource), itarget(ntarget)
        integer laddr(2,200)
        integer nterms(0:200)
        integer list(10 000)
        integer box(20)
        integer box1(20)
        integer itable(-3:3,-3:3,-3:3)
        integer nterms_eval(4,0:200)
        real *8 epsfmm,radius,t1,t2
        real *8 sourcesort(3,nsource)
        real *8 rlame(2)
        real *8 rlameuse(2)
        real *8 chargesort(3,nsource)
        real *8 dipstrsort(3,nsource)
        real *8 dipvecsort(3,nsource)
        real *8 rlams(nlams)
        real *8 whts(nlams)
        real *8 rlsc(nlams,0:nmax,0:nmax)
        real *8 skels(3,nskel,nboxes)
        real *8 targetsort(3,ntarget)
        real *8 wlists(lwlists)
        real *8 rmlexp(1)
        real *8 vec3(3,nchwork)
        real *8 quadvec(6,nchwork)
        real *8 timeinfo(10)
        real *8 center(3)
        real *8 scale(0:200)
        real *8 bsize(0:200)
        real *8 center0(3),corners0(3,8)
        real *8 center1(3),corners1(3,8)
        complex *16 ima
        complex *16 mexpf(nexptot)
        complex *16 mexpphys(nexptotp)
        complex *16 ptfrc(3)
        complex *16 hess(3,3)
        complex *16 fexpe(nexpe)
        complex *16 fexpback(nexpback)
        complex *16 chargeskel(nskel,nboxes)
        complex *16 ptfrctarg(3,ntarget)
        complex *16 hesstarg(3,3,ntarget)
        complex *16 chwork(nchwork)
        complex *16 mptemp(lmptemp)
        complex *16 mptemp2(lmptemp)
        complex *16 ptemp,ftmp(3),hesstmp(3,3)
        complex *16 ftmp2(3),hesstmp2(3,3)
c
        real *8 second,omp_get_wtime
c
        data ima/(0.0d0,1.0d0)/

c
c
c     INPUT PARAMETERS:
c
c     iprec        precision flag (see above)   ELIMINATE???
c     ifevalfar    far field flag (1 means compute far field, 
c                                  else dont)
c     ifevalloc    local field flag (1 means compute local field, 
c                                    else dont)
c     nsource      number of sources
c     sourcesort   sorted source coordinates
c     isource      sorting index for sources
c     ifsingle     flag indicating potential includes contribution
c                  from charges
c     chargesort   sorted charge values
c     ifdouble     flag indicating potential includes contribution
c                  from dipoles
c     dipstrsort   sorted dipole strengths
c     dipvecsort   sorted dipole orientation vectors
c     ifpot        potential flag (1 => compute, else do not)
c     iffld        field flag (1 => compute, else do not)
c     ntarget      number of targets
c     targetsort   sorted array of target locations
c     itarget      sorting index for targets
c     ifpottarg    target potential flag (1 => compute, else do not)
c     iffldtarg    target field flag (1 => compute, else do not)
c     epsfmm       FMM tolerance
c     iaddr        iaddr(2,nboxes) array points to mpole/local
c                     expansions for each box
c     rmlexp       workspace to contain mpole/local expansions.
c     xnodes       workspace to hold quadrature nodes
c     wts          workspace to hold quadrature weights
c     nquad        number of quadrature nodes
c     nboxes       number of boxes in FMM hierarchy
c     laddr        indexing array for FMM data structure
c     nlev         number of levels in FMM hierarchy
c     scale        array of scaling parameters
c     bsize        box dimension for FMM
c     nterms       array of nterms needed at each level
c     wlists       FMM data structure (real array)
c     lw           length of wlists
c
c
c     OUTPUT PARAMETERS:
c
c     pot          surface potential (if ifpot=1)
c     fld          surface field=-gradient(potential) (if iffld=1)
c     pottarg      target potential (if ifpot=1)
c     fldtarg      target field=-gradient(potential) (if iffld=1)
c     ier          error return code
c                  ier = 0    =>   normal execution
c                  ier = 4    =>   cannot allocate tree workspace
c                  ier = 8    =>   cannot alocate bulk FMM workspace
c                  ier = 16   =>   cannot allocate mpole exp workspace
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
      ifprint=1
c
c     
c       ... set the potential and field to zero
c
      do i=1,ntarget
         if( ifptfrctarg .eq. 1) then
            ptfrctarg(1,i)=0
            ptfrctarg(2,i)=0
            ptfrctarg(3,i)=0
         endif
         if( ifhesstarg .eq. 1) then
            do j = 1,3
            do k = 1,3
               hesstarg(j,k,i)=0
            enddo
            enddo
         endif
      enddo
c
      do i=1,10
         timeinfo(i)=0
      enddo
c
c
      if( ifevalfar .eq. 0 ) goto 8000
c       
c
c       ... initialize Legendre function evaluation routines
c
      nthmax = 0
      nthmaxp = 0
      do i=1,nlams
         if (numfour(i).gt.nthmax) nthmax = numfour(i)
         if (numphys(i).gt.nthmaxp) nthmaxp = numphys(i)
      enddo
c
      do i=0,nlev
      do itype=1,4
         call l3dterms_eval(itype,epsfmm,
     1        nterms_eval(itype,i),ier)
      enddo
      enddo
c
      if (ifprint .ge. 2) then
         call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
      endif
c
c       ... set all multipole and local expansions to zero
c
      do ibox = 1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
         level=box(1)
         call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
         call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
         call l3dzero(rmlexp(iaddr(3,ibox)),nterms(level))
         do j = 1,nexptotp*2
            rmlexp(iaddr(4,ibox)+j-1) = 0.0d0
            rmlexp(iaddr(5,ibox)+j-1) = 0.0d0
            rmlexp(iaddr(6,ibox)+j-1) = 0.0d0
         enddo
      enddo
c
      if (ifprint .ge. 1) call prinf('=== STEP 1 (form mp) ====*',i,0)
      t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
      do ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftmp,mptemp,chwork,quadvec,vec3) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            call d3tnkids(box,nkids)
c
            level=box(1)
c
c
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,20)
               call prinf('nkids=*',nkids,1)
            endif
c
            if (nkids .eq. 0) then
               npts=box(15)
               if (ifprint .ge. 2) then
                  call prinf('npts=*',npts,1)
                  call prinf('isource=*',isource(box(14)),box(15))
               endif
            endif
c
c       ... prune all sourceless boxes
c
            if( box(15) .eq. 0 ) cycle
c
            if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	       radius = (corners0(1,1) - center0(1))**2
	       radius = radius + (corners0(2,1) - center0(2))**2
	       radius = radius + (corners0(3,1) - center0(3))**2
	       radius = sqrt(radius)
c
               call formmp_mindlinb(scale(level),rlame,
     $          sourcesort(1,box(14)),ifsingle,chargesort(1,box(14)),
     $          ifdouble,dipstrsort(1,box(14)),dipvecsort(1,box(14)),
     $          npts,center0,nterms(level),mptemp,chwork,quadvec,
     $         rmlexp(iaddr(1,ibox)))
c
               call formmp_mindlinc(scale(level),rlame,
     $          sourcesort(1,box(14)),ifsingle,chargesort(1,box(14)),
     $          ifdouble,dipstrsort(1,box(14)),dipvecsort(1,box(14)),
     $          npts,center0,nterms(level),mptemp,chwork,quadvec,vec3,
     $          rmlexp(iaddr(2,ibox)),rmlexp(iaddr(3,ibox)))
c
            endif
c
         enddo
C$OMP END PARALLEL DO
      enddo
c
      t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
      timeinfo(1)=t2-t1
c       
      if (ifprint .ge. 1) call prinf('=== STEP 2 (form lo) ====*',i,0)
      t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and fields directly
c 
      do ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
         if( nlist .ne. 0 ) then
ccc         call prinf('nlist3=*', nlist,1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(level,npts,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftmp,ilist,mptemp,mptemp2) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
            do ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c  check no kids
c
               level1=box1(1)
               ifdirect3 = 0
ccc                  if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
ccc     $             box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
               if (box1(17).eq.0) cycle
c
               if( ifdirect3 .eq. 0 ) then
                  npts=box(15)
c
	          radius = (corners1(1,1) - center1(1))**2
	          radius = radius + (corners1(2,1) - center1(2))**2
	          radius = radius + (corners1(3,1) - center1(3))**2
	          radius = sqrt(radius)
c
c
                  call l3dformtaintkerbf90(ier,scale(level1),
     1               sourcesort(1,box(14)),ifsingle,
     1               chargesort(1,box(14)),ifdouble,
     $               dipstrsort(1,box(14)),dipvecsort(1,box(14)),
     $               npts,rlame,center1,radius,nterms(level1),mptemp)
                     call l3dadd(mptemp,rmlexp(iaddr(1,jbox)),
     $                    nterms(level1))
c
                  call l3dformtaintkercf90(ier,scale(level1),
     1               sourcesort(1,box(14)),ifsingle,
     1               chargesort(1,box(14)),ifdouble,
     $               dipstrsort(1,box(14)),dipvecsort(1,box(14)),
     $               npts,rlame,center1,radius,nterms(level1),mptemp,
     $               mptemp2)
                  call l3dadd(mptemp,rmlexp(iaddr(2,jbox)),
     $                    nterms(level1))
                  call l3dadd(mptemp2,rmlexp(iaddr(3,jbox)),
     $                    nterms(level1))
               else
                  call lfmm3dmindlinpart_direct(box,box1,rlame,
     $                  sourcesort,ifsingle,chargesort,ifdouble,
     $                  dipstrsort,dipvecsort,
     $                  targetsort,ifptfrctarg,ptfrctarg,
     $                  ifhesstarg,hesstarg)
               endif
            enddo
C$OMP END PARALLEL DO
         endif
      enddo
c
      t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
      timeinfo(2)=t2-t1
c
c
      if (ifprint .ge. 1) call prinf('=== STEPS 3,4,5 ====*',i,0)
      ifdirect2 = 0
      if (ifdirect2.eq.0) then
         call lfmm3d_list2_mindlin
     $        (bsize,rlame,nlev,laddr,scale,nterms,rmlexp,
     $        iaddr,nboxes,timeinfo,wlists,mptemp,mptemp2,lmptemp,
     $        nlams,numfour,nexptot,rlsc,nmax,mexpf,
     $        nthmax,fexpe,nexpe,fexpback,nexpback,rlams,whts,numphys,
     $        nthmaxp,nexptotp,skels,chargeskel,nskel)
      else
         do ibox=1,nboxes
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            itype=2
            call d3tgetl(ier,ibox,itype,list,nlist,wlists)
            if (nlist .gt. 0) then 
               if (ifprint .ge. 2) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('list2=*',list,nlist)
               endif
            endif
            if( box(15) .eq. 0 ) nlist=0
            if( nlist .ne. 0 ) then
               do ilist=1,nlist
                  jbox=list(ilist)
                  call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
                  call lfmm3dmindlinpart_direct(box,box1,rlame,
     $                  sourcesort,ifsingle,chargesort,ifdouble,
     $                  dipstrsort,dipvecsort,
     $                  targetsort,ifptfrctarg,ptfrctarg,
     $                  ifhesstarg,hesstarg)
               enddo
            endif
         enddo
      endif
c
c
      if (ifprint .ge. 1) call prinf('=== STEP 6 (eval mp) ====*',i,0)
      t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and fields directly
c
      nops=nops+1
      do ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         level = box(1)
c
c        compute equivalent charges
c
	 radius = (corners0(1,1) - center0(1))**2
	 radius = radius + (corners0(2,1) - center0(2))**2
	 radius = radius + (corners0(3,1) - center0(3))**2
	 radius = sqrt(radius)
         nquad = 2*nterms(level)
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if( box(15) .eq. 0 ) nlist=0
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
         if( nlist .ne. 0 ) then
ccc         call prinf('nlist4=*', nlist,1)
C$OMP PARALLEL DO reduction(+:nops) DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,level)
C$OMP$PRIVATE(ier,i,j,ptemp,ftmp,hesstmp,jj,ll,kk,ilist) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
            do ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
               level=box(1)
c
               ifdirect4 = 0
c
ccc               if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
ccc     $            box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
               if (ifdirect4 .eq. 0) then

ccc                  call prinf(' about to mindlinmpeval ibox is *',ibox,1)
ccc                  if (box1(17).gt.0) 
                  do jj = box1(16),box1(16)+box1(17)-1
                     call mindlinmpevalb(rlame,scale(level),center0,
     $                    rmlexp(iaddr(1,ibox)),nterms(level),
     $                    targetsort(1,jj),skels(1,1,ibox),
     $                    chargeskel(1,ibox),nskel,
     $                    ftmp,ifhesstarg,hesstmp)
                     if (ifptfrctarg .eq. 1) then
                        ptfrctarg(1,jj) = ptfrctarg(1,jj) +ftmp(1)
                        ptfrctarg(2,jj) = ptfrctarg(2,jj) +ftmp(2)
                        ptfrctarg(3,jj) = ptfrctarg(3,jj) +ftmp(3)
                     endif
                     if (ifhesstarg.eq.1) then
                        do kk = 1,3
                        do ll = 1,3
                           hesstarg(kk,ll,jj) = hesstarg(kk,ll,jj)+
     $                             hesstmp(kk,ll)
                        enddo
                        enddo
                     endif
                     call mindlinmpevalc(rlame,scale(level),center0,
     $                    rmlexp(iaddr(2,ibox)),rmlexp(iaddr(3,ibox)),
     $                    nterms(level),targetsort(1,jj),
     $                    ftmp,ifhesstarg,hesstmp)
                     if (ifptfrctarg.eq.1) then
                        ptfrctarg(1,jj) = ptfrctarg(1,jj)+
     $                         targetsort(3,jj)*ftmp(1)
                        ptfrctarg(2,jj) = ptfrctarg(2,jj)+
     $                         targetsort(3,jj)*ftmp(2)
                        ptfrctarg(3,jj) = ptfrctarg(3,jj)+
     $                         targetsort(3,jj)*ftmp(3)
                     endif
                     if (ifhesstarg.eq.1) then
                        do ll = 1,3
                           do kk = 1,3
                              hesstarg(kk,ll,jj)=hesstarg(kk,ll,jj)+
     $                           targetsort(3,jj)*hesstmp(kk,ll)
                           enddo
                           hesstarg(3,ll,jj) = hesstarg(3,ll,jj)-
     $                          ftmp(ll)
                        enddo
                     endif
                  enddo
               else
                  call lfmm3dmindlinpart_direct(box,box1,rlame,
     $                 sourcesort,ifsingle,chargesort,ifdouble,
     $                 dipstrsort,dipvecsort,
     $                 targetsort,ifptfrctarg,ptfrctarg,
     $                 ifhesstarg,hesstarg)
               endif
            enddo
C$OMP END PARALLEL DO
         endif
      enddo
c
ccc        call prinf('nops=*',nops,1)
c
      t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
      timeinfo(6)=t2-t1
c

      if (ifprint .ge. 1) call prinf('=== STEP 7 (eval lo) ====*',i,0)
      t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$PRIVATE(ftmp,hesstmp,jj,kk,ll)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
      do ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
         call d3tnkids(box,nkids)
c
ccc         ifprint = 0
         if (ifprint .ge. 2) then
            call prinf('ibox=*',ibox,1)
            call prinf('box=*',box,20)
            call prinf('nkids=*',nkids,1)
         endif
c
         if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
               call prinf('itarget=*',itarget(box(16)),box(17))
            endif
         endif
c
         if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
            level=box(1)
            npts=box(15)
ccc        call prinf('before taeval loop npts=*',npts,1)
c       
            if (level .ge. 2) then
               do jj = box(16),box(16)+box(17)-1
                  call mindlintaevalb(scale(level),center0,
     $                 rmlexp(iaddr(1,ibox)),nterms(level),
     $                 targetsort(1,jj),ftmp,ifhesstarg,hesstmp)
                  if (ifptfrctarg.eq.1) then
                     ptfrctarg(1,jj) = ptfrctarg(1,jj) +ftmp(1)
                     ptfrctarg(2,jj) = ptfrctarg(2,jj) +ftmp(2)
                     ptfrctarg(3,jj) = ptfrctarg(3,jj) +ftmp(3)
                  endif
                  if (ifhesstarg.eq.1) then
                     do kk = 1,3
                     do ll = 1,3
                        hesstarg(kk,ll,jj) = hesstarg(kk,ll,jj)+
     $                    hesstmp(kk,ll)
                     enddo
                     enddo
                  endif
                  call mindlintaevalc(scale(level),center0,
     $                 rmlexp(iaddr(2,ibox)),rmlexp(iaddr(3,ibox)),
     $                 nterms(level),targetsort(1,jj),
     $                 ftmp,ifhesstarg,hesstmp)
                  if (ifptfrctarg.eq.1) then
                     ptfrctarg(1,jj) = ptfrctarg(1,jj)+
     $                  targetsort(3,jj)*ftmp(1)
                     ptfrctarg(2,jj) = ptfrctarg(2,jj)+
     $                   targetsort(3,jj)*ftmp(2)
                     ptfrctarg(3,jj) = ptfrctarg(3,jj)+
     $                  targetsort(3,jj)*ftmp(3)
                  endif
                  if (ifhesstarg.eq.1) then
                     do ll = 1,3
                        do kk = 1,3
                           hesstarg(kk,ll,jj) = hesstarg(kk,ll,jj)+
     $                       targetsort(3,jj)*hesstmp(kk,ll)
                        enddo
                        hesstarg(3,ll,jj) = hesstarg(3,ll,jj)-
     $                        ftmp(ll)
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
C$OMP END PARALLEL DO
      t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
      timeinfo(7)=t2-t1
c
c
8000  continue
c
c
      if( ifevalloc .eq. 0 ) goto 9000
c 
      if (ifprint .ge. 1) call prinf('=== STEP 8 (direct) =====*',i,0)
      t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
      do ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
         call d3tnkids(box,nkids)
c
ccc        ifprint=2
         if (ifprint .ge. 2) then
            call prinf('ibox=*',ibox,1)
            call prinf('box=*',box,20)
            call prinf('nkids=*',nkids,1)
         endif
c
         if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
         endif
c
c
         if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
ccc        call prinf('calling directself *',nkids,1)
         call lfmm3dmindlinpart_direct_self(box,rlame,sourcesort,
     $        ifsingle,chargesort,ifdouble,dipstrsort,dipvecsort,
     $        targetsort,ifptfrctarg,ptfrctarg,ifhesstarg,hesstarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
            itype=1
            call d3tgetl(ier,ibox,itype,list,nlist,wlists)
            if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c
            do ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
               if( box1(15) .eq. 0 ) cycle
c    
               call lfmm3dmindlinpart_direct(box1,box,rlame,sourcesort,
     $            ifsingle,chargesort,ifdouble,dipstrsort,dipvecsort,
     $            targetsort,ifptfrctarg,ptfrctarg,ifhesstarg,hesstarg)
            enddo
         endif
      enddo
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
9000  continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
      ifprint=1
      if( ifprint .ge. 1 ) then
        call prinf
     $   ('formmp     formta     lolo      lota     tata     mpeval*',
     $   i,0)
        call prinf('taeval     direct*',i,0)
      endif
      if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
c       
      if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
      if (ifprint .ge. 1) call prinf('nsource=*',nsource,1)
      if (ifprint .ge. 1) call prinf('ntarget=*',ntarget,1)
c       
      return
      end
c
c
c
      subroutine lfmm3dmindlinpart_direct_self(box,rlame,
     $     source,ifsingle,charge,ifdouble,dipstr,dipvec,
     $     target,ifptfrctarg,ptfrctarg,ifhesstarg,hesstarg)
      implicit none
      integer ifsingle,ifdouble
      integer icomp,ifptfrctarg,ifhesstarg,j,i,ione
      integer box(20),box1(20)
c
      real *8 source(3,1),dipvec(3,1),rlame(2)
      real *8 charge(3,1),dipstr(3,1)
      real *8 target(3,1)
c
      complex *16 ptfrc(3,1),ptfrctarg(3,1)
      complex *16 hess(3,3,1),hesstarg(3,3,1)
      complex *16 ptemp,ftmp(3)
c
      if( ifptfrctarg .eq. 1 .or. ifhesstarg .eq. 1 ) then
         do j=box(16),box(16)+box(17)-1
         do i=box(14),box(14)+box(15)-1
            ione = 1
            do icomp = 1,3
               call intkerbc(icomp,rlame,
     $             source(1,i),ifsingle,charge(1,i),
     $             ifdouble,dipstr(1,i),dipvec(1,i),ione,
     $             target(1,j),ptemp,ifhesstarg,ftmp)
               if (ifptfrctarg.eq.1) then 
                  ptfrctarg(icomp,j) = ptfrctarg(icomp,j) + ptemp
               endif
               if (ifhesstarg.eq.1) then 
                  hesstarg(1,icomp,j) = hesstarg(1,icomp,j)+ftmp(1)
                  hesstarg(2,icomp,j) = hesstarg(2,icomp,j)+ftmp(2)
                  hesstarg(3,icomp,j) = hesstarg(3,icomp,j)+ftmp(3)
               endif
            enddo
         enddo
         enddo
      endif
      return
      end
c
c
c
c
      subroutine lfmm3dmindlinpart_direct(box,box1,rlame,
     $     source,ifsingle,charge,ifdouble,dipstr,dipvec,
     $     target,ifptfrctarg,ptfrctarg,ifhesstarg,hesstarg)
      implicit none
      integer ifsingle,ifdouble
      integer icomp,ifptfrctarg,ifhesstarg,j
      integer box(20),box1(20)
c
      real *8 source(3,*),dipvec(3,*),rlame(2)
      real *8 charge(3,*),dipstr(3,*)
      real *8 target(3,*)
c
      complex *16 ptfrctarg(3,*),hesstarg(3,3,*)
      complex *16 ptemp,ftmp(3)
c
      if( ifptfrctarg .eq. 1 .or. ifhesstarg .eq. 1 ) then
         do j=box1(16),box1(16)+box1(17)-1
            do icomp = 1,3
               call intkerbc(icomp,rlame,
     1            source(1,box(14)),ifsingle,charge(1,box(14)),
     2            ifdouble,dipstr(1,box(14)),dipvec(1,box(14)),box(15),
     2            target(1,j),ptemp,ifhesstarg,ftmp)
               if (ifptfrctarg.eq.1) then 
                  ptfrctarg(icomp,j) = ptfrctarg(icomp,j) + ptemp
               endif
               if (ifhesstarg.eq.1) then 
                  hesstarg(1,icomp,j) = hesstarg(1,icomp,j)+ftmp(1)
                  hesstarg(2,icomp,j) = hesstarg(2,icomp,j)+ftmp(2)
                  hesstarg(3,icomp,j) = hesstarg(3,icomp,j)+ftmp(3)
               endif
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c
c
      subroutine intkerbc(icomp,rlame,source,
     1           ifsingle,charge,ifdouble,dipstr,dipvec,
     2           nsource,ztrg,pot,iffld,fld)
      implicit none
      integer icomp,ifsingle,ifdouble,nsource,iffld
      integer i
      real *8 source(3,nsource),ztrg(3),charge(3,nsource)
      real *8 dipstr(3,nsource),dipvec(3,nsource),rlame(2)
      complex *16 cd,pot,fld(3)
      complex *16 potloc1,fldloc1(3)
      complex *16 potloc2,fldloc2(3)
c
c     computes field at target from a set of sources.
c
      pot = 0.0d0
      fld(1) = 0.0d0
      fld(2) = 0.0d0
      fld(3) = 0.0d0
      do i = 1,nsource
         call intker_mindlinb(icomp,rlame,source(1,i),ifsingle,
     1        charge(1,i),ifdouble,dipstr(1,i),dipvec(1,i),ztrg,
     1        potloc1,iffld,fldloc1)
         call intker_mindlinc(icomp,rlame,source(1,i),ifsingle,
     1        charge(1,i),ifdouble,dipstr(1,i),dipvec(1,i),ztrg,
     1        potloc2,iffld,fldloc2)
         pot = pot + potloc1 + ztrg(3)*potloc2
	 if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc1(1) + ztrg(3)*fldloc2(1)
         fld(2) = fld(2) + fldloc1(2) + ztrg(3)*fldloc2(2)
         fld(3) = fld(3) + fldloc1(3) + ztrg(3)*fldloc2(3)
         fld(3) = fld(3) - potloc2
         endif
      enddo        
      return
      end
c
c
c
c
c
c
c
c
      subroutine l3dpartmindlinalloc(wlists,iaddr,nboxes,lmptot,
     1            nterms,nexptotp)
      implicit none
      integer nboxes,lmptot,nexptotp
      integer iptr,ibox,ier,level
      integer box(20)
      integer iaddr(6,nboxes)
      integer nterms(0:*)
      integer wlists(*)
      real *8 center0(3),corners0(3,8)
c
c       ... construct pointer array iaddr for addressing multipole and
c       local expansion
c
      iptr=1
      do ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
         level=box(1)
c
c       ... first, allocate memory for three multipole expansions
c       
         iaddr(1,ibox)=iptr
         iptr=iptr+(nterms(level)+1)*(2*nterms(level)+1)*2
         iaddr(2,ibox)=iptr
         iptr=iptr+(nterms(level)+1)*(2*nterms(level)+1)*2
         iaddr(3,ibox)=iptr
         iptr=iptr+(nterms(level)+1)*(2*nterms(level)+1)*2
c
c       now add memory for three plane wave expansions
c
         iaddr(4,ibox)=iptr
         iptr=iptr+nexptotp*2
         iaddr(5,ibox)=iptr
         iptr=iptr+nexptotp*2
         iaddr(6,ibox)=iptr
         iptr=iptr+nexptotp*2
      enddo
      lmptot = iptr
      return
      end
c
c
c
c
c
c
c
      subroutine l3dpartreorderv2(nsource,source,
     $     ifsingle,charge,isource,ifdouble,
     1     dipstr,dipvec,sourcesort,chargesort,dipvecsort,dipstrsort) 
c
c     resort sources,etc. according to box location
c
      implicit none
      integer nsource,isource(nsource),ifsingle,ifdouble,i
      real *8 source(3,nsource),sourcesort(3,nsource)
      real *8 dipvec(3,nsource),dipvecsort(3,nsource)
      real *8 charge(3,nsource),chargesort(3,nsource)
      real *8 dipstr(3,nsource),dipstrsort(3,nsource)
c       
      do i = 1,nsource
         sourcesort(1,i) = source(1,isource(i))
         sourcesort(2,i) = source(2,isource(i))
         sourcesort(3,i) = source(3,isource(i))
         if( ifsingle .eq. 1 ) then
            chargesort(1,i) = charge(1,isource(i))
            chargesort(2,i) = charge(2,isource(i))
            chargesort(3,i) = charge(3,isource(i))
         endif
         if (ifdouble .eq. 1) then
            dipstrsort(1,i) = dipstr(1,isource(i))
            dipstrsort(2,i) = dipstr(2,isource(i))
            dipstrsort(3,i) = dipstr(3,isource(i))
            dipvecsort(1,i) = dipvec(1,isource(i))
            dipvecsort(2,i) = dipvec(2,isource(i))
            dipvecsort(3,i) = dipvec(3,isource(i))
         endif
      enddo
      return
      end
c
c
