cc Copyright (C) 2011: Leslie Greengard and Zydrunas Gimbutas
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
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the routines for Laplace particle FMM in R^3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lfmm3d_list2_mindlin
     $   (bsize,rlame,nlev,laddr,scale,nterms,rmlexp,iaddr,
     $   nboxes,timeinfo,wlists,mptemp,mptemp2,lmptemp,
     $   nlams,numfour,nexptot,rlsc,nmax,mexpf,
     $   nthmax,fexpe,nexpe,fexpback,nexpback,rlams,whts,numphys,
     $   nthmaxp,nexptotp,skels,chargeskel,nskel)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     List 2 processing routine for Mindlin B,C images
c
c     INPUT:
c
c     bsize      box sizes in FMM hierarchy
c     rlame      Lame coefficients
c     nlev       number of levels in FMM hierarchy
c                (0 to nlev) or (1 to nlev+1) in integer
c                arrays beginning at 1 rather than 0.
c     laddr      indexing array for FMM data structure
c     scale      array of scaling parameters
c     nterms     array of nterms needed at each level
c     rmlexp     array of outgoing mpole/planewave expansions.
c     rmlexp2    array of incoming mpole/planewave expansions.
c     iaddr      indexing array for mpole/planewave expansions
c     nboxes     number of boxes in FMM data structure
c     timeinfo   timing data IN/OUT
c     wlists     FMM data structure (real array)
c     mptemp     workspace for mpole expansion
c     lmptemp    length of mptemp
C     nlams      number of discret. pts. in plane wave lambda integral 
C     numfour    number of Fourier modes for each lambda
C     nexptot    sum_j numfour(j)
C     rlsc       precomputed array needed by mtx_mindlin routines
C     nmax       dimension parameter for rlsc
C     mexpf      workspace for Fourier rep plane wave expansion
C     nthmax      max_j numfour(j)
C     fexpe      precomputed array needed by mtx_mindlin routines
C     nexpe      length of fexpe
C     fexpback   precomputed array needed by phystof routines
C     nexpback   length of fexpback
C     rlams      discretization points in lambda integral 
C     whts       discretization weights in lambda integral 
C     numphys    number of physical space modes for each lambda
C     nthmaxp    max_j numphys(j)
C     nexptotp   sum_j numphys(j)
C     nskel      number of pts in discret of charge ring (see paper).
C
C     OUTPUT:
C
C     skels      charge ring discretization pts for each box.
C     chargeskel charge ring strengths     
C     rmlexp     multipole B expansions are overwritten by 
C                 diff between original expansions and expansions
C                 due to skels (very complicated analysis - see paper).
C     rmlexp2     contains all shifted expansion data 
C-----------------------------------------------------------------------
ccc      implicit real *8 (a-h,o-z)
      implicit none
c
      integer iaddr(6,nboxes),laddr(2,nlev+1),nterms(0:nlev)
      integer nlev,lmptemp,nlams,nexptot,nthmax,nexpe,nthmaxp,nexptotp
      integer nboxes,nmax,nexpback,nskel,i,ibox,ier,ifdirect2
      integer ifprint,ilev,ilist,itype,jbox,level,level0,level1
      integer nkids,nlist,nquad
      integer numfour(nlams)
      integer numphys(nlams)
      integer list(10000)
      integer box(20)
      integer box1(20)
      real *8 rlame(2),radius,scc,t1,t2,xdis,ydis,zdis
      real *8 rlams(nlams),whts(nlams)
      real *8 rlsc(nlams,0:nmax,0:nmax)
      real *8 skels(3,nskel,nboxes)
      real *8 scale(0:nlev)
      real *8 rmlexp(*)
      real *8 bsize(0:nlev)
      real *8 center0(3),corners0(3,8)
      real *8 center1(3),corners1(3,8)
      real *8 wlists(*)
      complex *16 chargeskel(nskel,nboxes)
      complex *16 fexpe(nexpe)
      complex *16 fexpback(nexpe)
      complex *16 mexpf(nexptot)
      complex *16 mptemp(lmptemp)
      complex *16 mptemp2(lmptemp)
      complex *16 ptemp,ftemp(3)
      real *8, allocatable :: rlampow(:)
      real *8, allocatable :: facts(:)
      complex *16, allocatable :: zeye(:)
c      real *8 rlampow(0:nmax)
c      real *8 facts(0:2*nmax)
c      complex *16 zeye(0:nthmax)
c
      real *8 timeinfo(10)
      real *8 second, omp_get_wtime
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
      ifprint=0
c      allocate(rlampow(0:nmax))
c      allocate(facts(0:2*nmax))
c      allocate(zeye(0:nthmax))
c
c
      if (ifprint .ge. 1) then
         call prinf('=== STEP 3 (merge mp) ===*',i,0)
      endif
      t1=second()
C$    t1=omp_get_wtime()
c
c       ... step 3, merge all multipole expansions
c       
      do ilev=nlev,3,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(ier,i,ptemp,ftemp,scc) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            call d3tnkids(box,nkids)
c
c       ... prune all sourceless boxes
c
            if ((box(15) .gt. 0) .and. (nkids .ne.0) ) then
c
               level0=box(1)
               if( level0 .ge. 2 ) then
                  radius = (corners0(1,1) - center0(1))**2
                  radius = radius + (corners0(2,1) - center0(2))**2
                  radius = radius + (corners0(3,1) - center0(3))**2
                  radius = sqrt(radius)
c       
                  if( ifprint .ge. 2 ) then
                     call prin2('radius=*',radius,1)
                     call prinf('ibox=*',ibox,1)
                     call prinf('box=*',box,20)
                     call prinf('nkids=*',nkids,1)
                  endif
c
c       ... merge multipole expansions of the kids 
c
                  if (ifprint .ge. 2) then
                      call prin2('center0=*',center0,3)
                  endif
                  do i = 1,8
                     jbox = box(5+i)
                     if (jbox.gt.0) then
                       call d3tgetb
     $                 (ier,jbox,box1,center1,corners1,wlists)
                       if (ifprint .ge. 2) then
                          call prinf('jbox=*',jbox,1)
                          call prin2('center1=*',center1,3)
                       endif
                       level1=box1(1)
                       if (box1(15).ne.0) then
                         call l3dmpmpquadu_add(scale(level1),center1,
     $                    rmlexp(iaddr(1,jbox)),nterms(level1),
     $                    scale(level0),center0,rmlexp(iaddr(1,ibox)),
     $                    nterms(level0),nterms(level0),ier)
                         call l3dmpmpquadu_add(scale(level1),center1,
     $                    rmlexp(iaddr(2,jbox)),nterms(level1),
     $                    scale(level0),center0,rmlexp(iaddr(2,ibox)),
     $                    nterms(level0),nterms(level0),ier)
                         call l3dmpmpquadu_add(scale(level1),center1,
     $                    rmlexp(iaddr(3,jbox)),nterms(level1),
     $                    scale(level0),center0,rmlexp(iaddr(3,ibox)),
     $                    nterms(level0),nterms(level0),ier)
                       endif
                     endif
                  enddo
c
c     We are using the Laplace eq mpmp shift. 
c     The rescaling is off by a constant factor for Mindlin B.
c     mprstemp corrects this.
c
                  scc = scale(level1)/scale(level0)
                  call  mprstemp(rmlexp(iaddr(1,ibox)),
     $                  nterms(level0),scc)
                  if (ifprint .ge. 2) then
                     call prinf('=============*',scc,0)
                  endif
               endif
            endif
         enddo
C$OMP END PARALLEL DO
      enddo
c
        call prinf('=== UPWARD PASS COMPLETE ===*',i,0)
c
c
      t2=second()
C$    t2=omp_get_wtime()
      timeinfo(3)=t2-t1
c
      if (ifprint .ge. 1) then
          call prinf('=== STEP 4 (mp to lo) ===*',i,0)
      endif
      t1=second()
C$    t1=omp_get_wtime()
c
c       ... step 4, convert multipole expansions into the local ones
c
      do ilev=3,nlev+1

c       t3=second()
cC$      t3=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0)
C$OMP$PRIVATE(ier,nquad,radius)
C$OMP$PRIVATE(mptemp,mexpf)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)  
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,20)
            endif
            level0=box(1)
            if (box(15).gt.0) then
               call mtx_mindlinb(scale(level0),rlame,
     $           rmlexp(iaddr(1,ibox)),nterms(level0),nlams,numfour,
     $           nexptot,mexpf,rlsc,nthmax,fexpe,nexpe,rlams,whts,
     $           numphys,nthmaxp,nexptotp,rmlexp(iaddr(4,ibox)))
c
               nquad = 2*nmax
               radius = (corners0(1,1) - center0(1))**2
               radius = radius + (corners0(2,1) - center0(2))**2
               radius = radius + (corners0(3,1) - center0(3))**2
               radius = sqrt(radius)
               call mptoslpborder(rmlexp(iaddr(1,ibox)),nterms(level0),
     $           center0,scale(level0),radius,nquad,skels(1,1,ibox),
     $           chargeskel(1,ibox))
               call l3dformmp(ier,scale(level0),skels(1,1,ibox),
     $           chargeskel(1,ibox),nskel,center0,nterms(level0),
     $           mptemp)
               call mpdiff(mptemp,rmlexp(iaddr(1,ibox)),nterms(level0))
               call mpdzminus2(mptemp,nterms(level0))
               call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level0))
               call l3dadd(mptemp,rmlexp(iaddr(1,ibox)),
     $           nterms(level0))

               call mtx_mindlinc(rlame,rmlexp(iaddr(2,ibox)),
     $           rmlexp(iaddr(3,ibox)),nterms(level0),nlams,numfour,
     $           nexptot,mexpf,rlsc,nthmax,fexpe,nexpe,rlams,whts,
     $           numphys,nthmaxp,nexptotp,rmlexp(iaddr(5,ibox)),
     $           rmlexp(iaddr(6,ibox)))
            endif
         enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(ier,i,ilist,mptemp,mptemp2,xdis,ydis,zdis)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)  
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,20)
            endif
            level0=box(1)
            if (level0 .ge. 2) then
c
c       ... retrieve list #2
c
               itype=2
               call d3tgetl(ier,ibox,itype,list,nlist,wlists)
               if (ifprint .ge. 2) then
                  call prinf('list2=*',list,nlist)
               endif
c
c       ... for all pairs in list #2, apply the translation operator 
c
               do ilist=1,nlist
                  jbox=list(ilist)
                  call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
                  if ((jbox.gt.0) .and. 
     $             (box1(15).gt.0) .and. 
     $             (box(17).gt.0)) then
                     radius = (corners1(1,1) - center1(1))**2
                     radius = radius + (corners1(2,1) - center1(2))**2
                     radius = radius + (corners1(3,1) - center1(3))**2
                     radius = sqrt(radius)
c
c       ... convert multipole expansions for all boxes in list 2 in local exp
c       ... if source is childless, evaluate directly (if cheaper)
c
                     level1=box1(1)
                     ifdirect2 = 0
                     zdis = (center0(3)-center1(3))/bsize(level1)
                     if (abs(zdis).lt. 1.5) then
                        ifdirect2 = 1
                     endif
c
                     if (ifdirect2 .eq. 0) then
                       xdis = (center0(1)-center1(1))/bsize(level1)
                       ydis = (center0(2)-center1(2))/bsize(level1)
                       call xyzshift_add(rmlexp(iaddr(4,jbox)),
     $                   rlams,nlams,numphys,rmlexp(iaddr(4,ibox)),
     $                   nexptotp,xdis,ydis,zdis)
                       call xyzshift_add(rmlexp(iaddr(5,jbox)),
     $                   rlams,nlams,numphys,rmlexp(iaddr(5,ibox)),
     $                   nexptotp,xdis,ydis,zdis)
                       call xyzshift_add(rmlexp(iaddr(6,jbox)),
     $                   rlams,nlams,numphys,rmlexp(iaddr(6,ibox)),
     $                   nexptotp,xdis,ydis,zdis)
                     else
                       call mindlinbmploc(scale(level0),center1,
     $    	         rmlexp(iaddr(1,jbox)),nterms(level0),radius,
     $                   scale(level1),center0,mptemp,
     $                   nterms(level1),skels(1,1,jbox),
     $                   chargeskel(1,jbox),nskel,mptemp2,rlame)
                       call l3dadd(mptemp,rmlexp(iaddr(1,ibox)),
     $                   nterms(level1))
	               call l3dmplocquadu_add(scale(level0),center1,
     $    	         rmlexp(iaddr(2,jbox)),nterms(level0),
     $                   scale(level1),center0,rmlexp(iaddr(2,ibox)),
     $                   nterms(level1),nterms(level1),ier)
	               call l3dmplocquadu_add(scale(level0),center1,
     $    	         rmlexp(iaddr(3,jbox)),nterms(level0),
     $                   scale(level1),center0,rmlexp(iaddr(3,ibox)),
     $                   nterms(level1),nterms(level1),ier)
                     endif
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO

c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0)
C$OMP$PRIVATE(ier,nquad,radius)
C$OMP$PRIVATE(mptemp,mexpf)
C$OMP$PRIVATE(rlampow,facts,zeye)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)  
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,20)
            endif
            level0=box(1)
            if (box(17).gt.0) then
      allocate(rlampow(0:nmax))
      allocate(facts(0:2*nmax))
      allocate(zeye(0:nthmax))
               call phystof(mexpf,nexptot,nlams,rlams,numfour,numphys,
     $            nthmax,rmlexp(iaddr(4,ibox)),nexptotp,fexpback,
     $            nexpback)
               call exptolocald_add(rmlexp(iaddr(1,ibox)),nterms,
     $            rlams,whts,nlams,numfour,nthmax,nexptot,mexpf,mptemp,
     $            bsize(level0),scale(level0),rlampow,facts,zeye)
               call phystof(mexpf,nexptot,nlams,rlams,numfour,numphys,
     $             nthmax,rmlexp(iaddr(5,ibox)),nexptotp,fexpback,
     $             nexpback)
               call exptolocald_add(rmlexp(iaddr(2,ibox)),nterms,
     $            rlams,whts,nlams,numfour,nthmax,nexptot,mexpf,mptemp,
     $            bsize(level0),scale(level0),rlampow,facts,zeye)
               call phystof(mexpf,nexptot,nlams,rlams,numfour,numphys,
     $             nthmax,rmlexp(iaddr(6,ibox)),nexptotp,fexpback,
     $             nexpback)
               call exptolocald_add(rmlexp(iaddr(3,ibox)),nterms,
     $            rlams,whts,nlams,numfour,nthmax,nexptot,mexpf,mptemp,
     $            bsize(level0),scale(level0),rlampow,facts,zeye)
      deallocate(rlampow)
      deallocate(facts)
      deallocate(zeye)
            endif
         enddo
C$OMP END PARALLEL DO

c        t4=second()
cC$      t4=omp_get_wtime()

      enddo
c
      t2=second()
C$    t2=omp_get_wtime()
      timeinfo(4)=t2-t1
c       
      if (ifprint .ge. 1) then
         call prinf('=== STEP 5 (split lo) ===*',i,0)
      endif
      t1=second()
C$    t1=omp_get_wtime()
c
c       ... step 5, split all local expansions
c
      do ilev=3,nlev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(ier,i) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
            call d3tgetb(ier,ibox,box,center0,corners0,wlists)
            call d3tnkids(box,nkids)
c       
            if (nkids .ne. 0) then
               level0=box(1)
               if (level0 .ge. 2) then
                  if (ifprint .ge. 2) then
                     call prinf('ibox=*',ibox,1)
                     call prinf('box=*',box,20)
                     call prinf('nkids=*',nkids,1)
                     call prin2('center0=*',center0,3)
                  endif
c
c       ... split local expansion of the parent box
c
                  do i = 1,8
	            jbox = box(5+i)
ccc	            if (jbox.gt.0) then
	            if ((jbox.gt.0) .and. (box(17).gt.0)) then
                      call d3tgetb(ier,jbox,box1,center1,corners1,
     $                     wlists)
                      radius = (corners1(1,1) - center1(1))**2
                      radius = radius + (corners1(2,1) - center1(2))**2
                      radius = radius + (corners1(3,1) - center1(3))**2
                      radius = sqrt(radius)
                      if (ifprint .ge. 2) then
                         call prinf('jbox=*',jbox,1)
                         call prin2('radius=*',radius,1)
                         call prin2('center1=*',center1,3)
                      endif
                      level1=box1(1)
	              call l3dloclocquadu_add(scale(level0),center0,
     $    	         rmlexp(iaddr(1,ibox)),nterms(level0),
     $                   scale(level1),center1,rmlexp(iaddr(1,jbox)),
     $                   nterms(level1),nterms(level1),ier)
	              call l3dloclocquadu_add(scale(level0),center0,
     $    	         rmlexp(iaddr(2,ibox)),nterms(level0),
     $                   scale(level1),center1,rmlexp(iaddr(2,jbox)),
     $                   nterms(level1),nterms(level1),ier)
	              call l3dloclocquadu_add(scale(level0),center0,
     $    	         rmlexp(iaddr(3,ibox)),nterms(level0),
     $                   scale(level1),center1,rmlexp(iaddr(3,jbox)),
     $                   nterms(level1),nterms(level1),ier)
                    endif
                  enddo
                  if (ifprint .ge. 2) then 
                     call prinf('=============*',scc,0)
                  endif
               endif
            endif
c
            if (nkids .ne. 0) then
                level=box(1)
                if (level .ge. 2) then
                   if( ifprint .ge. 2 ) then
                      call prinf('ibox=*',ibox,1)
                      call prinf('box=*',box,20)
                      call prinf('nkids=*',nkids,1)
                   endif
                endif
            endif
         enddo
C$OMP END PARALLEL DO
      enddo
c       
      t2=second()
C$        t2=omp_get_wtime()
      timeinfo(5)=t2-t1
      return
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mprstemp(mpole,nterms,sc)
c
c     rescale multipole expansion by factor sc.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i,j,nterms
      real *8 sc
      complex *16 mpole(0:nterms,-nterms:nterms)
c
      do i = 0,nterms
      do j = -nterms,nterms
         mpole(i,j) = mpole(i,j)*sc
      enddo
      enddo
      return
      end
