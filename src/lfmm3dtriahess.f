cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-07-10 19:11:13 -0400 (Sat, 10 Jul 2010) $
c    $Revision: 1079 $
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Laplace potentials, fields, and
c     hessians on surfaces defined by a collection of positively
c     oriented flat triangles.  (FORTRAN 90 VERSION)
c
c     lfmm3dtriahess - Laplace FMM in R^3: evaluate all pairwise triangle
c         interactions (including self-interaction)
c
c     lfmm3dtriahesstarg - Laplace FMM in R^3: evaluate all pairwise
c         triangle interactions (including self-interaction) +
c         interactions with targets
c
c     l3dtriahessdirecttarg - Laplace interactions in R^3: evaluate all pairwise
c         triangle interactions (including self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
c     The routines in this file permit the calculation of SECOND 
c     DERIVATIVES of the potentials. ARBITRARY ORIENTED dipole vectors 
c     are permitted.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine lfmm3dtriahess(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the harmonic potential, field, and
c       hessian due to a collection of flat triangles with constant
c       single and/or double layer densities. We use (1/r) for the
c       Green's function, without the (1/4 pi) scaling.
c       Self-interactions are included.
c   
c       The main FMM routine permits both evaluation on surface
c       and at a collection of off-surface targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c
        lused7=0
c
        ntarget=0
        ifpottarg=0
        iffldtarg=0
        ifhesstarg=0
c
        call lfmm3dtriahesstarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriahesstarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       Laplace FMM in R^3: evaluate all pairwise triangle
c       interactions (including self-interaction) + interactions with targets
c
c       This is the principal subroutine for evaluating 
c       harmonic layer potentials on (flat) triangulated surfaces.
c       It permits the evaluation of a single layer potential
c       with piecewise constant density defined by <<charge>>
c       and a dipole layer potential with piecewise constant density 
c       and dipole orientation defined by <<dipstr,dipvec>>.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are included.
c   
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dtriahesstargmain.
c
c       NOTE: In this routine, arbitrary oriented dipole vectors are permitted.
c
c
c       INPUT PARAMETERS:
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
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): piecewise constant SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): piecewise constant DLP strengths
c       dipvec: real *8 (3,nsource): piecewise constant dipole orientation 
c                                    vectors. 
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
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
c       pot: complex *16 (nsource): potential at triangle centroids
c       fld: complex *16 (3,nsource): field (-gradient) at triangle centroids 
c       hess: complex *16 (6,nsource): hessian at triangle centroids 
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c       hesstarg: complex *16 (6,ntarget): hessian at target locations
c
cf2py   intent(out) ier
cf2py   intent(in) iprec
cf2py   intent(in) nsource
cf2py   intent(in) triaflat,trianorm,source
cf2py   intent(in) ifcharge,charge
cf2py   check(!ifcharge || (shape(charge,0) == nsource))  charge
cf2py   depend(nsource)  charge
cf2py   intent(in) ifdipole,dipvec,dipstr
cf2py   check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
cf2py   depend(nsource)  dipstr
cf2py   intent(in) ifpot,iffld,ifhess
cf2py   intent(out) pot,fld,hess
cf2py   intent(in) ifpottarg, iffldtarg, ifhesstarg
cf2py   intent(in) target
cf2py   intent(in) ntarget
cf2py   check((!ifpottarg && !iffldtarg && !ifhesstarg) || (shape(target,0)==3 && shape(target,1) == ntarget))  target
cf2py   check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
cf2py   check((!iffldtarg) || (shape(fldtarg,0)==3 && shape(fldtarg,1) == ntarget))  fldtarg
cf2py   check((!ifhesstarg) || (shape(hesstarg,0)==6 && shape(hesstarg,1) == ntarget))  hesstarg
c
c       (F2PY workaround: *targ must be input because f2py
c       refuses to allocate zero-size output arrays.)
c
cf2py   intent(in,out) pottarg,fldtarg, hesstarg

        dimension triaflat(3,3,nsource)
        dimension trianorm(3,nsource)
        dimension source(3,nsource)
        dimension dipvec(3,nsource)
        dimension target(3,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        complex *16 ima
        complex *16 pot(nsource)
        complex *16 fld(3,nsource)
        complex *16 hess(6,nsource)
        complex *16 pottarg(ntarget)
        complex *16 fldtarg(3,ntarget)
        complex *16 hesstarg(6,ntarget)
c
        dimension timeinfo(10)
c       
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 levels of refinement, which is 
c     more than enough for any non-pathological case.
c
c       
        dimension laddr(2,200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:) 
        real *8, allocatable :: wlists(:) 
        real *8, allocatable :: wrmlexp(:) 
        complex *16 ptemp,ftemp(3)
c       
        data ima/(0.0d0,1.0d0)/
c              
        ier=0
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
        ifprint=1
c
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
        if (ifprint.eq.1) call prin2('epsfmm=*',epsfmm,1)
c
c     set criterion for box subdivision (number of sources per box)
c
        if( iprec .eq. -2 ) nbox=8/3
        if( iprec .eq. -1 ) nbox=15/3
        if( iprec .eq. 0 ) nbox=30/3
        if( iprec .eq. 1 ) nbox=60/3
        if( iprec .eq. 2 ) nbox=120/3
        if( iprec .eq. 3 ) nbox=240/3
        if( iprec .eq. 4 ) nbox=480/3
        if( iprec .eq. 5 ) nbox=700/3
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        if (ifprint.eq.1) call prinf('nbox=*',nbox,1)
c
c     create oct-tree data structure
c
        ntot = 100*(nsource+ntarget)+10000
        if( iprec .eq. -2 ) ntot = ntot * 1.5*1.5*1.5
        if( iprec .eq. -1 ) ntot = ntot * 1.5*1.5
        do ii = 1,10
           allocate (wlists(ntot))
           call lfmm3dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
              goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return          
        endif
c
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
c       
        if (ifprint.eq.1) call prin2('scale=*',scale,nlev+1)
c
c       carve up workspace further
c
c     itriaflatsort is pointer for sorted triangle coordinates
c     itrianormsort is pointer for sorted triangle normals
c     isourcesort is pointer for sorted triangle centroids
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
c
        itriaflatsort = lused7 + 5
        ltriaflatsort = 3*3*nsource
        itrianormsort = itriaflatsort + ltriaflatsort
        ltrianormsort = 3*nsource
        isourcesort = itrianormsort + ltrianormsort 
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 3*nsource
          ldipstr = 2*nsource
        else
          ldipvec = 3
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr       
c
c     allocate workspace for potential and field arrays
c
        ipot = lused7
        lpot = 2*nsource
        lused7=lused7+lpot
c       
        ifld = lused7
        if( iffld .eq. 1) then
           lfld = 2*(3*nsource)
        else
           lfld=6
        endif
        lused7=lused7+lfld
c      
        ihess = lused7
        if( ifhess .eq. 1) then
        lhess = 2*(6*nsource)
        else
        lhess=12
        endif
        lused7=lused7+lhess
c      
        ipottarg = lused7
        lpottarg = 2*ntarget
        lused7=lused7+lpottarg
c       
        ifldtarg = lused7
        if( iffldtarg .eq. 1) then
           lfldtarg = 2*(3*ntarget)
        else
           lfldtarg=6
        endif
        lused7=lused7+lfldtarg
c      
        ihesstarg = lused7
        if( ifhesstarg .eq. 1) then
        lhesstarg = 2*(6*ntarget)
        else
        lhesstarg=12
        endif
        lused7=lused7+lhesstarg
c      
        if (ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c      
        nmax = 0
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call l3dterms(epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        nquad=2*nmax        
c       
c     ixnodes is pointer for quadrature nodes
c     iwhts is pointer for quadrature weights
c
        ixnodes = lused7 
        iwts = ixnodes + nquad
        lused7 = iwts + nquad
c
        if (ifprint.eq.1) call prinf('nterms=*',nterms,nlev+1)
        if (ifprint.eq.1) call prinf('nmax=*',nmax,1)
c
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c   
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7 
        imptemp = iiaddr + 2*nboxes
        lmptemp = (nmax+1)*(2*nmax+1)*2 
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                   lused7 is *',lused7,1)
           ier = 8
           return          
        endif
c
c     reorder triangles, centroids etc. so that each box holds
c     contiguous list of source numbers.
c
        call l3dreorder(nsource,source,ifcharge,charge,wlists(iisource),
     $     ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 
c       
        call l3dreordertria(nsource,wlists(iisource),
     $     triaflat,w(itriaflatsort),trianorm,w(itrianormsort))
c
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
      if (ifprint.eq.1) then
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
      endif
c       
        ifinit=1
        call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
c
ccc        call prin2('xnodes=*',xnodes,nquad)
ccc        call prin2('wts=*',wts,nquad)
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call l3dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if (ifprint.eq.1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                   lused7 is *',lused7,1)
           ier = 16
           return          
        endif
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
c
c
c
        call lfmm3dtriahesstargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,w(itriaflatsort),w(itrianormsort),
     $     w(isourcesort),w(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),iffld,w(ifld),ifhess,w(ihess),
     $     ntarget,w(itargetsort),w(iitarget),
     $     ifpottarg,w(ipottarg),iffldtarg,w(ifldtarg),
     $     ifhesstarg,w(ihesstarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     w(ixnodes),w(iwts),nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if(ifprint.eq.1) then
        call prinf('lwlists=*',lwlists,1)
        call prinf('lused total =*',lused7,1)       
        call prin2('memory / point = *',(lused7)/dble(nsource),1)
        endif
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
        if(iffld .eq. 1) 
     $     call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
        if(ifhess .eq. 1) 
     $     call l3dhsort(nsource,wlists(iisource),w(ihess),hess)
c
        if(ifpottarg .eq. 1 )
     $     call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
        if(iffldtarg .eq. 1) 
     $     call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
        if(ifhesstarg .eq. 1) 
     $     call l3dhsort(ntarget,wlists(iitarget),w(ihesstarg),hesstarg)
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriahesstargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,triaflatsort,trianormsort,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ifhess,hess,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,xnodes,wts,nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension triaflatsort(3,3,1),trianormsort(3,1)
        dimension sourcesort(3,1), isource(1)
        complex *16 chargesort(1)
        complex *16 dipstrsort(1)
        dimension dipvecsort(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1)
        dimension targetsort(3,1), itarget(1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        complex *16 hesstarg(6,1)
        dimension wlists(1)
        dimension iaddr(2,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension xnodes(nquad),wts(nquad)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3),htemp(6)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(40 000)
        dimension nterms_eval(4,0:200)
c
        real *8, allocatable :: scarray_local(:)
        real *8, allocatable :: scarray_mpole(:)
c
        data ima/(0.0d0,1.0d0)/

c
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
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(1,i)=0
           hess(2,i)=0
           hess(3,i)=0
           hess(4,i)=0
           hess(5,i)=0
           hess(6,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(1,i)=0
           hesstarg(2,i)=0
           hesstarg(3,i)=0
           hesstarg(4,i)=0
           hesstarg(5,i)=0
           hesstarg(6,i)=0
        endif
        enddo
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        norder=1
        nqtri=1
c
        if( iprec .eq. -2 ) then
        norder=2
        nqtri=2
        endif
c
        if( iprec .eq. -1 ) then
        norder=2
        nqtri=2
        endif
c
        if( iprec .eq. 0 ) then
        norder=4
        nqtri=4
        endif
c
        if( iprec .ge. 1 ) then
        norder=6
        nqtri=6
        endif
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c       ... initialize Legendre function evaluation routines
c
        nlege=100
        lw7=40 000
        call ylgndrfwini(nlege,wlege,lw7,lused7)
c
        do i=0,nlev
        do itype=1,4
        call l3dterms_eval(itype,epsfmm,
     1       nterms_eval(itype,i),ier)
        enddo
        enddo
c
        if(ifprint .ge. 2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
        do ibox = 1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        enddo
c
c
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
        do 1300 ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
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
c        ipts=box(14)
c        npts=box(15)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(15)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(14)),box(15))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(15) .eq. 0 ) goto 1200
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
 	    if (ifcharge .eq. 1) then
               call l3dformmptris2_add(ier,scale(level),
     1  	  triaflatsort(1,1,box(14)),chargesort(box(14)),npts,
     2            center0,norder,nterms(level),rmlexp(iaddr(1,ibox)),
     $            wlege,nlege)
            endif
c
            if (ifdipole .eq. 1 ) then
               call l3dformmptrid2_hess_add(ier,scale(level),
     1           triaflatsort(1,1,box(14)),trianormsort(1,box(14)),
     2           dipstrsort(box(14)),dipvecsort(1,box(14)),
     $           npts,center0,norder,
     3           nterms(level),rmlexp(iaddr(1,ibox)),wlege,nlege)
            endif

         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
        if(ifprint .ge. 1) 
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and fields directly
c 
         do 3251 ibox=1,nboxes
c
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
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            level1=box1(1)
c
            ifdirect3 = 0
            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
c
               npts=box(15)
c
               if( ifcharge .eq. 1 ) then
                  call l3dformtatris2_add(ier,scale(level1),
     1  	    triaflatsort(1,1,box(14)),chargesort(box(14)),
     $              npts,center1,
     2              norder,nterms(level1),
     $              rmlexp(iaddr(2,jbox)),wlege,nlege)
	       endif
               if( ifdipole .eq. 1 ) then
                  call l3dformtatrid2_hess_add(ier,scale(level1),
     1               triaflatsort(1,1,box(14)),trianormsort(1,box(14)),
     2               dipstrsort(box(14)),dipvecsort(1,box(14)),
     $               npts,center1,norder,nterms(level1),
     $               rmlexp(iaddr(2,jbox)),wlege,nlege)
               endif
c
            else

            call lfmm3dtriahess_direct(nqtri,box,box1,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
            endif
         enddo
C$OMP END PARALLEL DO
c
 3251    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c
        if(ifprint .ge. 1) 
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0
        if (ifhess.eq.1) ifprune_list2 = 0
        call lfmm3d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2)
c
c
        allocate( scarray_mpole(0:100000) )
        call l3dmpevalhessdini(nterms(0),scarray_mpole)
c       
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and fields directly
c
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,level,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box(1)
c
            ifdirect4 = 0
c
            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then
               do j=box1(14),box1(14)+box1(15)-1
                  if( ifhess .eq. 1 ) then
c                  call l3dmpevalhess(scale(level),center0,
c     $               rmlexp(iaddr(1,ibox)),nterms(level),
c     $               sourcesort(1,j),
c     $               ptemp,iffld,ftemp,ifhess,htemp,
c     $               ier)
                  call l3dmpevalhessd_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),
     $               sourcesort(1,j),
     $               ptemp,iffld,ftemp,ifhess,htemp,
     $               scarray_mpole,wlege,nlege)
                  else
                  call l3dmpeval_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),nterms(level),
     $               sourcesort(1,j),
     $               ptemp,iffld,ftemp,
     $               wlege,nlege,ier)
                  endif
                  if( ifpot .eq. 1 ) pot(j)=pot(j)+ptemp
                  if( iffld .eq. 1 ) then
                     fld(1,j)=fld(1,j)+ftemp(1)
                     fld(2,j)=fld(2,j)+ftemp(2)
                     fld(3,j)=fld(3,j)+ftemp(3)
                  endif
                  if (ifhess .eq. 1) then
                  hess(1,j)=hess(1,j)+htemp(1)
                  hess(2,j)=hess(2,j)+htemp(2)
                  hess(3,j)=hess(3,j)+htemp(3)
                  hess(4,j)=hess(4,j)+htemp(4)
                  hess(5,j)=hess(5,j)+htemp(5)
                  hess(6,j)=hess(6,j)+htemp(6)
                  endif
               enddo
               do j=box1(16),box1(16)+box1(17)-1
                  if( ifhesstarg .eq. 1 ) then
c                  call l3dmpevalhess(scale(level),center0,
c     $               rmlexp(iaddr(1,ibox)),nterms(level),
c     $               targetsort(1,j),
c     $               ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
c     $               ier)
                  call l3dmpevalhessd_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),
     $               targetsort(1,j),
     $            ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
     $            scarray_mpole,wlege,nlege)
                  else
                  call l3dmpeval_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),nterms(level),
     $               targetsort(1,j),
     $               ptemp,iffldtarg,ftemp,
     $               wlege,nlege,ier)
                  endif
                  if( ifpottarg .eq. 1 ) pottarg(j)=pottarg(j)+ptemp
                  if( iffldtarg .eq. 1 ) then
                     fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                     fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                     fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
                  endif
                  if (ifhesstarg .eq. 1) then
                  hesstarg(1,j)=hesstarg(1,j)+htemp(1)
                  hesstarg(2,j)=hesstarg(2,j)+htemp(2)
                  hesstarg(3,j)=hesstarg(3,j)+htemp(3)
                  hesstarg(4,j)=hesstarg(4,j)+htemp(4)
                  hesstarg(5,j)=hesstarg(5,j)+htemp(5)
                  hesstarg(6,j)=hesstarg(6,j)+htemp(6)
                  endif
               enddo
            else
            
            call lfmm3dtriahess_direct(nqtri,box,box1,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
            endif
        enddo
C$OMP END PARALLEL DO
 3252   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c
        allocate( scarray_local(0:100000) )
        call l3dtaevalhessdini(nterms(0),scarray_local)

        if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
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
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(15)
c       
        if (level .ge. 2) then
c
            do j=box(14),box(14)+box(15)-1
               if( ifhess .eq. 1 ) then
c                  call l3dtaevalhess(scale(level),center0,
c     $               rmlexp(iaddr(2,ibox)),nterms(level),
c     $               sourcesort(1,j),
c     $               ptemp,iffld,ftemp,ifhess,htemp,
c     $               ier)
               call l3dtaevalhessd_trunc(scale(level),center0,
     $              rmlexp(iaddr(2,ibox)),nterms(level),sourcesort(1,j),
     $              ptemp,iffld,ftemp,ifhess,htemp,
     $              scarray_local,wlege,nlege)
               else
               call l3dtaeval_trunc(scale(level),center0,
     $            rmlexp(iaddr(2,ibox)),nterms(level),nterms(level),
     $            sourcesort(1,j),
     $            ptemp,iffld,ftemp,
     $            wlege,nlege,ier)
               endif
               if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
               if (iffld .eq. 1) then
                  fld(1,j)=fld(1,j)+ftemp(1)
                  fld(2,j)=fld(2,j)+ftemp(2)
                  fld(3,j)=fld(3,j)+ftemp(3)
               endif
               if (ifhess .eq. 1) then
               hess(1,j)=hess(1,j)+htemp(1)
               hess(2,j)=hess(2,j)+htemp(2)
               hess(3,j)=hess(3,j)+htemp(3)
               hess(4,j)=hess(4,j)+htemp(4)
               hess(5,j)=hess(5,j)+htemp(5)
               hess(6,j)=hess(6,j)+htemp(6)
               endif
            enddo

            do j=box(16),box(16)+box(17)-1
               if( ifhesstarg .eq. 1 ) then
c                  call l3dtaevalhess(scale(level),center0,
c     $               rmlexp(iaddr(2,ibox)),nterms(level),
c     $               targetsort(1,j),
c     $               ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
c     $               ier)
               call l3dtaevalhessd_trunc(scale(level),center0,
     $              rmlexp(iaddr(2,ibox)),nterms(level),
     $              targetsort(1,j),
     $              ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
     $              scarray_local,wlege,nlege)
               else
               call l3dtaeval_trunc(scale(level),center0,
     $            rmlexp(iaddr(2,ibox)),nterms(level),nterms(level),
     $            targetsort(1,j),
     $            ptemp,iffldtarg,ftemp,
     $            wlege,nlege,ier)
               endif
               if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
               if (iffldtarg .eq. 1) then
                  fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                  fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                  fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
               endif
               if (ifhesstarg .eq. 1) then
               hesstarg(1,j)=hesstarg(1,j)+htemp(1)
               hesstarg(2,j)=hesstarg(2,j)+htemp(2)
               hesstarg(3,j)=hesstarg(3,j)+htemp(3)
               hesstarg(4,j)=hesstarg(4,j)+htemp(4)
               hesstarg(5,j)=hesstarg(5,j)+htemp(5)
               hesstarg(6,j)=hesstarg(6,j)+htemp(6)
               endif
            enddo
        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
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
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
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
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        call lfmm3dtriahess_direct_self(nqtri,box,
     $     triaflatsort,trianormsort,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
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
c    
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c
            call lfmm3dtriahess_direct(nqtri,box1,box,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
c
 6203       continue

        endif
c
 6202   continue
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
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if(ifprint .ge. 1) then
        call prin2('timeinfo=*',timeinfo,8)       
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        endif
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriahess_direct(nqtri,box,box1,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),hess(6,1)
        complex *16 pottarg(1),fldtarg(3,1),hesstarg(6,1)
        complex *16 ptemp,ftemp(3),htemp(6)
c
c
c       ... sources
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 .or. ifhess .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(14),box1(14)+box1(15)-1
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps_hess3(box(15),source(1,j),
     $     charge(box(14)),triaflat(1,1,box(14)),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
        hess(1,j)=hess(1,j)+htemp(1)
        hess(2,j)=hess(2,j)+htemp(2)
        hess(3,j)=hess(3,j)+htemp(3)
        hess(4,j)=hess(4,j)+htemp(4)
        hess(5,j)=hess(5,j)+htemp(5)
        hess(6,j)=hess(6,j)+htemp(6)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd_hess3(box(15),source(1,j),
     $     dipstr(box(14)),dipvec(1,box(14)),triaflat(1,1,box(14)),
     $     trianorm(1,box(14)),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
        hess(1,j)=hess(1,j)+htemp(1)
        hess(2,j)=hess(2,j)+htemp(2)
        hess(3,j)=hess(3,j)+htemp(3)
        hess(4,j)=hess(4,j)+htemp(4)
        hess(5,j)=hess(5,j)+htemp(5)
        hess(6,j)=hess(6,j)+htemp(6)
        endif
        endif
c
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 
     $     .or. ifhesstarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(16),box1(16)+box1(17)-1
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps_hess3(box(15),target(1,j),
     $     charge(box(14)),triaflat(1,1,box(14)),
     $     ifpottarg,ptemp,iffldtarg,ftemp,ifhesstarg,htemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
        hesstarg(1,j)=hesstarg(1,j)+htemp(1)
        hesstarg(2,j)=hesstarg(2,j)+htemp(2)
        hesstarg(3,j)=hesstarg(3,j)+htemp(3)
        hesstarg(4,j)=hesstarg(4,j)+htemp(4)
        hesstarg(5,j)=hesstarg(5,j)+htemp(5)
        hesstarg(6,j)=hesstarg(6,j)+htemp(6)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd_hess3(box(15),target(1,j),
     $     dipstr(box(14)),dipvec(1,box(14)),triaflat(1,1,box(14)),
     $     trianorm(1,box(14)),
     $     ifpottarg,ptemp,iffldtarg,ftemp,ifhesstarg,htemp)        
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
        hesstarg(1,j)=hesstarg(1,j)+htemp(1)
        hesstarg(2,j)=hesstarg(2,j)+htemp(2)
        hesstarg(3,j)=hesstarg(3,j)+htemp(3)
        hesstarg(4,j)=hesstarg(4,j)+htemp(4)
        hesstarg(5,j)=hesstarg(5,j)+htemp(5)
        hesstarg(6,j)=hesstarg(6,j)+htemp(6)
        endif
        endif
c
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriahess_direct_self(nqtri,box,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),hess(6,1)
        complex *16 pottarg(1),fldtarg(3,1),hesstarg(6,1)
        complex *16 ptemp,ftemp(3),htemp(6)
c
        ione = 1
c
c       ... sources
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 
     $     .or. ifhess .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(14),box(14)+box(15)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifcharge.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilaps_hess3(ione,ione,
     1     source(1,j),charge(j),triaflat(1,1,j),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        else
        call direct3dtritarglaps_hess3(ione,source(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
        hess(1,j)=hess(1,j)+htemp(1)
        hess(2,j)=hess(2,j)+htemp(2)
        hess(3,j)=hess(3,j)+htemp(3)
        hess(4,j)=hess(4,j)+htemp(4)
        hess(5,j)=hess(5,j)+htemp(5)
        hess(6,j)=hess(6,j)+htemp(6)
        endif
        endif
c
        if (ifdipole.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilapd_hess3(ione,ione,
     $     source(1,j),dipstr(j),dipvec(1,j),triaflat(1,1,j),
     $     trianorm(1,j),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        else
        call direct3dtritarglapd_hess3(ione,source(1,j),
     $     dipstr(i),dipvec(1,i),triaflat(1,1,i),
     $     trianorm(1,i),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
        hess(1,j)=hess(1,j)+htemp(1)
        hess(2,j)=hess(2,j)+htemp(2)
        hess(3,j)=hess(3,j)+htemp(3)
        hess(4,j)=hess(4,j)+htemp(4)
        hess(5,j)=hess(5,j)+htemp(5)
        hess(6,j)=hess(6,j)+htemp(6)
        endif
        endif
c
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 
     $     .or. ifhesstarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(16),box(16)+box(17)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps_hess3(ione,target(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp,ifhesstarg,htemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
        hesstarg(1,j)=hesstarg(1,j)+htemp(1)
        hesstarg(2,j)=hesstarg(2,j)+htemp(2)
        hesstarg(3,j)=hesstarg(3,j)+htemp(3)
        hesstarg(4,j)=hesstarg(4,j)+htemp(4)
        hesstarg(5,j)=hesstarg(5,j)+htemp(5)
        hesstarg(6,j)=hesstarg(6,j)+htemp(6)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd_hess3(ione,target(1,j),
     $     dipstr(i),dipvec(1,i),triaflat(1,1,i),
     $     trianorm(1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp,ifhesstarg,htemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
        hesstarg(1,j)=hesstarg(1,j)+htemp(1)
        hesstarg(2,j)=hesstarg(2,j)+htemp(2)
        hesstarg(3,j)=hesstarg(3,j)+htemp(3)
        hesstarg(4,j)=hesstarg(4,j)+htemp(4)
        hesstarg(5,j)=hesstarg(5,j)+htemp(5)
        hesstarg(6,j)=hesstarg(6,j)+htemp(6)
        endif
        endif
c
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
        subroutine l3dtriahessdirecttarg(nsource,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,ntarget,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       Laplace interactions in R^3: evaluate all pairwise triangle
c       interactions (including self-interaction) + interactions with targets
c       cia direct O(N^2) algorithm.
c
c       This is the principal subroutine for evaluating 
c       harmonic layer potentials on (flat) triangulated surfaces.
c       It permits the evaluation of a single layer potential
c       with piecewise constant density defined by <<charge>>
c       and a dipole layer potential with piecewise constant density 
c       and dipole orientation defined by <<dipstr,dipvec>>.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are included.
c   
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dtriahesstargmain.
c
c       NOTE: In this routine, arbitrary oriented dipole vectors are permitted.
c
c
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): piecewise constant SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): piecewise constant DLP strengths
c       dipvec: real *8 (3,nsource): piecewise constant dipole orientation 
c                                    vectors. 
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at triangle centroids
c       fld: complex *16 (3,nsource): field (-gradient) at triangle centroids 
c       hess: complex *16 (6,nsource): hessian at triangle centroids 
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c       hesstarg: complex *16 (6,ntarget): hessian at target locations
c
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),hess(6,1)
        complex *16 pottarg(1),fldtarg(3,1),hesstarg(6,1)
        complex *16 ptemp,ftemp(3),htemp(6)
c
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(1,i)=0
           hess(2,i)=0
           hess(3,i)=0
           hess(4,i)=0
           hess(5,i)=0
           hess(6,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(1,i)=0
           hesstarg(2,i)=0
           hesstarg(3,i)=0
           hesstarg(4,i)=0
           hesstarg(5,i)=0
           hesstarg(6,i)=0
        endif
        enddo
c
        ione = 1
c
c       ... sources
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 
     $     .or. ifhess .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
        do j=1,nsource
        do i=1,nsource
c
        if (ifcharge.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilaps_hess3(ione,ione,
     1     source(1,j),charge(j),triaflat(1,1,j),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        else
        call direct3dtritarglaps_hess3(ione,source(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
        hess(1,j)=hess(1,j)+htemp(1)
        hess(2,j)=hess(2,j)+htemp(2)
        hess(3,j)=hess(3,j)+htemp(3)
        hess(4,j)=hess(4,j)+htemp(4)
        hess(5,j)=hess(5,j)+htemp(5)
        hess(6,j)=hess(6,j)+htemp(6)
        endif
        endif
c
        if (ifdipole.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilapd_hess3(ione,ione,
     $     source(1,j),dipstr(j),dipvec(1,j),triaflat(1,1,j),
     $     trianorm(1,j),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        else
        call direct3dtritarglapd_hess3(ione,source(1,j),
     $     dipstr(i),dipvec(1,i),triaflat(1,1,i),
     $     trianorm(1,i),
     $     ifpot,ptemp,iffld,ftemp,ifhess,htemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
        hess(1,j)=hess(1,j)+htemp(1)
        hess(2,j)=hess(2,j)+htemp(2)
        hess(3,j)=hess(3,j)+htemp(3)
        hess(4,j)=hess(4,j)+htemp(4)
        hess(5,j)=hess(5,j)+htemp(5)
        hess(6,j)=hess(6,j)+htemp(6)
        endif
        endif
c
        enddo
        enddo
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 
     $     .or. ifhesstarg .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
        do j=1,ntarget
        do i=1,nsource
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps_hess3(ione,target(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp,ifhesstarg,htemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
        hesstarg(1,j)=hesstarg(1,j)+htemp(1)
        hesstarg(2,j)=hesstarg(2,j)+htemp(2)
        hesstarg(3,j)=hesstarg(3,j)+htemp(3)
        hesstarg(4,j)=hesstarg(4,j)+htemp(4)
        hesstarg(5,j)=hesstarg(5,j)+htemp(5)
        hesstarg(6,j)=hesstarg(6,j)+htemp(6)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd_hess3(ione,target(1,j),
     $     dipstr(i),dipvec(1,i),triaflat(1,1,i),
     $     trianorm(1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp,ifhesstarg,htemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
        hesstarg(1,j)=hesstarg(1,j)+htemp(1)
        hesstarg(2,j)=hesstarg(2,j)+htemp(2)
        hesstarg(3,j)=hesstarg(3,j)+htemp(3)
        hesstarg(4,j)=hesstarg(4,j)+htemp(4)
        hesstarg(5,j)=hesstarg(5,j)+htemp(5)
        hesstarg(6,j)=hesstarg(6,j)+htemp(6)
        endif
        endif
c
        enddo
        enddo
C$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c***********************************************************************
c
c                    Laplace SLP routines
c
c***********************************************************************
      subroutine direct3dtrilaps_hess3(ipatch,ntri,
     1     zparts,charge,triang,
     $     ifpot,rpot,iffld,ptfrc,ifhess,hess)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise-constant charge density on collection of triangles 
c     numbered jpatch = 1,...,ntri. 
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c     
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     charge(ntri)      array of (piecewise constant) SLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c     hess(6)           hessian of potential at centroid
c
c
      implicit none
      integer ifpot,iffld,ifhess
      integer ifinit, ipatch, jpatch, ntri, nquad, ier
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      real *8   valxx,valyy,valzz,valxy,valxz,valyz
      real *8   derxx,deryy,derzz,derxy,derxz,deryz
      complex *16 charge(1), rpot, ptfrc(3), hess(6)
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      hess(1) = 0.0d0
      hess(2) = 0.0d0
      hess(3) = 0.0d0
      hess(4) = 0.0d0
      hess(5) = 0.0d0
      hess(6) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         if (ipatch.eq.jpatch) then
            iquad = 0
            if( ifpot .eq. 1 ) then
            itype = 1
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            endif
            if( iffld .eq. 1 ) then
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            valz = -valz
            endif
            if( ifhess .eq. 1 ) then
            itype = 12
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 14
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 13
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            valxz=-valxz
            valyz=-valyz
            endif
         else
            iquad = 0
            if (z0.gt.0) iquad = 1
            if (z0.lt.0) iquad = -1
            if( ifpot .eq. 1 ) then
            itype = 1
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            endif
            if( iffld .eq. 1 ) then
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            valz = -valz
            endif
            if( ifhess .eq. 1 ) then
            itype = 12
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 14
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 13
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            valxz=-valxz
            valyz=-valyz
            endif
         endif
c
         if( ifpot .eq. 1 ) then
         rpot = rpot + charge(jpatch)*val
         endif
c
         if( iffld .eq. 1 ) then
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*derx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*dery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*derz
         endif
c
         if( ifhess .eq. 1 ) then
         call rothess3d(w,triang(1,1,jpatch),
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)
         hess(1)=hess(1)+charge(jpatch)*derxx
         hess(2)=hess(2)+charge(jpatch)*deryy
         hess(3)=hess(3)+charge(jpatch)*derzz
         hess(4)=hess(4)+charge(jpatch)*derxy
         hess(5)=hess(5)+charge(jpatch)*derxz
         hess(6)=hess(6)+charge(jpatch)*deryz
         endif
c
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarglaps_hess3(ntri,targ,charge,triang,
     1           ifpot,rpot,iffld,ptfrc,ifhess,hess)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG not
c     lying on the surface due to piecewise-constant charge density
c     on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     charge(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c     hess(6)           hessian of potential at targ
c
c
c
      implicit none
      integer ifpot,iffld,ifhess
      integer ifinit, jpatch, ntri, nquad, ier,j
      integer iquad, itype
      real *8   targ(3), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
      real *8   valxx,valyy,valzz,valxy,valxz,valyz
      real *8   derxx,deryy,derzz,derxy,derxz,deryz
      complex *16 charge(1), rpot, ptfrc(3), hess(6)
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      hess(1) = 0.0d0
      hess(2) = 0.0d0
      hess(3) = 0.0d0
      hess(4) = 0.0d0
      hess(5) = 0.0d0
      hess(6) = 0.0d0
c
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         iquad = 0
         if (z0.gt.0) iquad = 1
         if (z0.lt.0) iquad = -1

         if( ifpot .eq. 1 ) then
         itype = 1
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         endif

         if( iffld .eq. 1 ) then
         itype = 2
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valx)
         itype = 3
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valy)
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valz)
         valz = -valz
         endif
c
         if( ifhess .eq. 1 ) then
         itype = 12
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valxx)
         itype = 14
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valyy)
         itype = 7
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valzz)
         itype = 13
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valxy)
         itype = 5
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valxz)
         itype = 6
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valyz)
         valxz=-valxz
         valyz=-valyz
         endif
c
         if( ifpot .eq. 1 ) then
         rpot = rpot + charge(jpatch)*val
         endif
c
         if( iffld .eq. 1 ) then
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*derx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*dery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*derz
         endif
c
         if( ifhess .eq. 1 ) then
         call rothess3d(w,triang(1,1,jpatch),
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)
         hess(1)=hess(1)+charge(jpatch)*derxx
         hess(2)=hess(2)+charge(jpatch)*deryy
         hess(3)=hess(3)+charge(jpatch)*derzz
         hess(4)=hess(4)+charge(jpatch)*derxy
         hess(5)=hess(5)+charge(jpatch)*derxz
         hess(6)=hess(6)+charge(jpatch)*deryz
         endif
c
      enddo
      return
      end
c
c
c***********************************************************************
c
c                 Laplace DLP routines
c
c***********************************************************************
      subroutine direct3dtrilapd_hess(ipatch,ntri,
     1           zparts,dipstr,triang,trinorm,rpot,ptfrc,hess)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     dipstr(ntri)      array of (piecewise constant) DLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c     hess(6)           hessian of potential at centroid
c
c
      implicit none
      integer ifinit, ipatch, jpatch, ntri, nquad, ier
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      real *8   valxx,valyy,valzz,valxy,valxz,valyz
      real *8   derxx,deryy,derzz,derxy,derxz,deryz
      complex *16 dipstr(1), rpot, ptfrc(3), hess(6)
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      hess(1) = 0.0d0
      hess(2) = 0.0d0
      hess(3) = 0.0d0
      hess(4) = 0.0d0
      hess(5) = 0.0d0
      hess(6) = 0.0d0
c
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
         if (ipatch.eq.jpatch) then
            itype = 4
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            itype = 5
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)

            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 20
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)

         else
            iquad = 0
            if (z0.gt.0) iquad = +1
            if (z0.lt.0) iquad = -1
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)

            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 20
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)

         endif

         rpot = rpot + dipstr(jpatch)*val
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
c
         call rothess3d(w,triang(1,1,jpatch),
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)
         hess(1)=hess(1)+dipstr(jpatch)*derxx
         hess(2)=hess(2)+dipstr(jpatch)*deryy
         hess(3)=hess(3)+dipstr(jpatch)*derzz
         hess(4)=hess(4)+dipstr(jpatch)*derxy
         hess(5)=hess(5)+dipstr(jpatch)*derxz
         hess(6)=hess(6)+dipstr(jpatch)*deryz
c
      enddo 
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarglapd_hess
     $     (ntri,targ,dipstr,triang,trinorm,rpot,ptfrc,hess)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     dipstr(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c     hess(6)           hessian of potential at targ
c
c
      implicit none
      integer ifinit, jpatch, ntri, nquad, ier,j
      integer iquad, itype
      real *8   targ(3), triang(3,3,1), trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
      real *8   valxx,valyy,valzz,valxy,valxz,valyz
      real *8   derxx,deryy,derzz,derxy,derxz,deryz
      complex *16 dipstr(1), rpot, ptfrc(3), hess(6)
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      hess(1) = 0.0d0
      hess(2) = 0.0d0
      hess(3) = 0.0d0
      hess(4) = 0.0d0
      hess(5) = 0.0d0
      hess(6) = 0.0d0
c
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         val = -val
         itype = 5
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valx)
         itype = 6
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valy)
         itype = 7
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valz)

            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 20
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)

         rpot = rpot + dipstr(jpatch)*val
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
c
         call rothess3d(w,triang(1,1,jpatch),
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)
         hess(1)=hess(1)+dipstr(jpatch)*derxx
         hess(2)=hess(2)+dipstr(jpatch)*deryy
         hess(3)=hess(3)+dipstr(jpatch)*derzz
         hess(4)=hess(4)+dipstr(jpatch)*derxy
         hess(5)=hess(5)+dipstr(jpatch)*derxz
         hess(6)=hess(6)+dipstr(jpatch)*deryz
c       
      enddo
      return
      end
c
c
c
c
c
      subroutine direct3dtrilapd_hess3(ipatch,ntri,
     1     zparts,dipstr,dipvec,triang,trinorm,
     $     ifpot,rpot,iffld,ptfrc,ifhess,hess)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     dipstr(ntri)      array of (piecewise constant) DLP strengths
c     dipvec(3,ntri)      array of (piecewise constant) DLP vectors
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c     hess(6)           hessian of potential at centroid
c
c
      implicit none
      integer ifpot,iffld,ifhess
      integer ifinit, ipatch, jpatch, ntri, nquad, ier, k
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   trinorm(3,1), dipvec(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12), vectout(3)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      real *8   valxx,valyy,valzz,valxy,valxz,valyz
      real *8   derxx,deryy,derzz,derxy,derxz,deryz
      complex *16 dipstr(1), rpot, ptfrc(3), hess(6)
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      hess(1) = 0.0d0
      hess(2) = 0.0d0
      hess(3) = 0.0d0
      hess(4) = 0.0d0
      hess(5) = 0.0d0
      hess(6) = 0.0d0
c
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
ccc         call tri_for_vect(w,trinorm(1,jpatch),vectout)
ccc         call prin2('after rotation, trinorm=*',vectout,3)

         call tri_for_vect(w,dipvec(1,jpatch),vectout)
ccc         call prin2('after rotation, dipvec=*',vectout,3)

         do k=1,3

         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1
         if (ipatch.eq.jpatch) iquad = 0

           if ( k .eq. 1 ) then

            if( ifpot .eq. 1 ) then
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif

            if( iffld .eq. 1 ) then
            itype = 12
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 13
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif

            if( ifhess .eq. 1 ) then
            itype = 21
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 23
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 22
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            endif

            endif


           if ( k .eq. 2 ) then

            if( ifpot .eq. 1 ) then
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif

            if( iffld .eq. 1 ) then
            itype = 13
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 14
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif

            if( ifhess .eq. 1 ) then
            itype = 22
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 24
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 23
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            endif

            endif


            if ( k .eq. 3 ) then

            if( ifpot .eq. 1 ) then
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif

            if( iffld .eq. 1 ) then
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif

            if( ifhess .eq. 1 ) then
            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 20
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            endif

            endif

            
         if( ifpot .eq. 1 ) then
         val=val*vectout(k)
         rpot = rpot + dipstr(jpatch)*val
         endif
         
         if( iffld .eq. 1 ) then
         valx=valx*vectout(k)
         valy=valy*vectout(k)
         valz=valz*vectout(k)
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
         endif
c
         if( ifhess .eq. 1 ) then
         valxx=valxx*vectout(k)
         valyy=valyy*vectout(k)
         valzz=valzz*vectout(k)
         valxy=valxy*vectout(k)
         valxz=valxz*vectout(k)
         valyz=valyz*vectout(k)
         call rothess3d(w,triang(1,1,jpatch),
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)
         hess(1)=hess(1)+dipstr(jpatch)*derxx
         hess(2)=hess(2)+dipstr(jpatch)*deryy
         hess(3)=hess(3)+dipstr(jpatch)*derzz
         hess(4)=hess(4)+dipstr(jpatch)*derxy
         hess(5)=hess(5)+dipstr(jpatch)*derxz
         hess(6)=hess(6)+dipstr(jpatch)*deryz
         endif
c
         enddo
c
      enddo 
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarglapd_hess3
     $     (ntri,targ,dipstr,dipvec,triang,trinorm,
     $     ifpot,rpot,iffld,ptfrc,ifhess,hess)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     dipstr(ntri)      array of SLP strengths (constant)
c     dipvec(3,ntri)      array of (piecewise constant) DLP vectors
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c     hess(6)           hessian of potential at targ
c
c
      implicit none
      integer ifpot,iffld,ifhess
      integer ifinit, jpatch, ntri, nquad, ier,j, k
      integer iquad, itype
      real *8   targ(3), triang(3,3,1), trinorm(3,1), dipvec(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12), vectout(3)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
      real *8   valxx,valyy,valzz,valxy,valxz,valyz
      real *8   derxx,deryy,derzz,derxy,derxz,deryz
      complex *16 dipstr(1), rpot, ptfrc(3), hess(6)
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      hess(1) = 0.0d0
      hess(2) = 0.0d0
      hess(3) = 0.0d0
      hess(4) = 0.0d0
      hess(5) = 0.0d0
      hess(6) = 0.0d0
c
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
ccc         call tri_for_vect(w,trinorm(1,jpatch),vectout)
ccc         call prin2('after rotation, trinorm=*',vectout,3)

         call tri_for_vect(w,dipvec(1,jpatch),vectout)
ccc         call prin2('after rotation, dipvec=*',vectout,3)
c
         do k=1,3

         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1

           if ( k .eq. 1 ) then

            if( ifpot .eq. 1 ) then
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif

            if( iffld .eq. 1 ) then
            itype = 12
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 13
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif

            if( ifhess .eq. 1 ) then
            itype = 21
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 23
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 22
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            endif

            endif


           if ( k .eq. 2 ) then

            if( ifpot .eq. 1 ) then
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif

            if( iffld .eq. 1 ) then
            itype = 13
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 14
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif

            if( ifhess .eq. 1 ) then
            itype = 22
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 24
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 23
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            endif

            endif


            if ( k .eq. 3 ) then

            if( ifpot .eq. 1 ) then
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif

            if( iffld .eq. 1 ) then
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif

            if( ifhess .eq. 1 ) then
            itype = 15
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxx)
            itype = 17
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyy)
            itype = 20
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valzz)
            itype = 16
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxy)
            itype = 18
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valxz)
            itype = 19
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valyz)
            endif

            endif

            
         if( ifpot .eq. 1 ) then
         val=val*vectout(k)
         rpot = rpot + dipstr(jpatch)*val
         endif
         
         if( iffld .eq. 1 ) then
         valx=valx*vectout(k)
         valy=valy*vectout(k)
         valz=valz*vectout(k)
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
         endif
c
         if( ifhess .eq. 1 ) then
         valxx=valxx*vectout(k)
         valyy=valyy*vectout(k)
         valzz=valzz*vectout(k)
         valxy=valxy*vectout(k)
         valxz=valxz*vectout(k)
         valyz=valyz*vectout(k)
         call rothess3d(w,triang(1,1,jpatch),
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)
         hess(1)=hess(1)+dipstr(jpatch)*derxx
         hess(2)=hess(2)+dipstr(jpatch)*deryy
         hess(3)=hess(3)+dipstr(jpatch)*derzz
         hess(4)=hess(4)+dipstr(jpatch)*derxy
         hess(5)=hess(5)+dipstr(jpatch)*derxz
         hess(6)=hess(6)+dipstr(jpatch)*deryz
         endif
c
         enddo
c
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
c
c       Form a multipole expansion due to arbitrary oriented constant
c       dipole density on a flat triangle
c
c***********************************************************************
c
c
c
      subroutine l3dformmptridone2_hess
     $     (ier, scale, triang, trinorm,
     1          dipstr, dipvec, x0y0z0, norder, nterms, mpole,
     $     wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c     Form multipole expansion due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal
c     dipstr       = DLP strength
c     dipvec       = dipole vector
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c     work(lw)        = workspace
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), trinorm(3), x0y0z0(3), dipvec(3)
      real *8  scale, wlege(1)
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  dipvec_interp(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20)
c
      complex *16 eye, dipstr, zk, dipstr_interp(500)
      complex *16  mpole(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc        call prinf('nnodes=*',nnodes,1)
c
c        dipvec_interp(1,1)=trinorm(1)
c        dipvec_interp(2,1)=trinorm(2)
c        dipvec_interp(3,1)=trinorm(3)
c
        dipvec_interp(1,1)=dipvec(1)
        dipvec_interp(2,1)=dipvec(2)
        dipvec_interp(3,1)=dipvec(3)

ccc        write(*,*) trinorm, dipvec
        
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         dipstr_interp(i) = dipstr*weights(i)
         call tri_bak(w,vertout,zparts(1,i))
         dipvec_interp(1,i) = dipvec_interp(1,1)
         dipvec_interp(2,i) = dipvec_interp(2,1)
         dipvec_interp(3,i) = dipvec_interp(3,1)
      enddo
      call l3dformmp_dp_trunc(ier,scale,zparts,dipstr_interp,
     $   dipvec_interp,nnodes,x0y0z0,nterms,nterms,mpole,wlege,nlege)
      return
      end
c
c
c
c
c
      subroutine l3dformmptrid2_hess_add(ier,scale,triang,trinorm,
     $     dipstr,dipvec,
     1     ntri,x0y0z0,norder,nterms,mpole,wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     dipvec             = array of dipole vectors
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c     w(lw)              = additional work array
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale, wlege(1)
      real *8  dipvec(3,1)
      complex *16  eye, zk, dipstr(1)
      complex *16  mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformmptridone2_hess
     $   (ier,scale,triang(1,1,i),trinorm(1,i),dipstr(i),dipvec(1,i),
     $   x0y0z0, norder,nterms,mptemp,wlege,nlege)
         call l3dadd(mptemp,mpole,nterms)
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
c
c       Form a local expansion due to arbitrary oriented constant
c       dipole density on a flat triangle
c
c***********************************************************************
      subroutine l3dformtatridone2_hess
     $     (ier, scale, triang, trinorm,
     1          dipstr, dipvec, x0y0z0, norder, nterms, local,
     $     wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c     Form local expansion due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal
c     dipstr       = DLP strength
c     dipvec       = dipole vector
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of local expansion
c     work(lw)        = workspace
c
c     OUTPUT:
c
c     local        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), trinorm(3), x0y0z0(3), dipvec(3)
      real *8  scale, wlege(1)
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  dipvec_interp(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20)
c
      complex *16 eye, dipstr, zk, dipstr_interp(500)
      complex *16  local(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc        call prinf('nnodes=*',nnodes,1)
c
c        dipvec_interp(1,1)=trinorm(1)
c        dipvec_interp(2,1)=trinorm(2)
c        dipvec_interp(3,1)=trinorm(3)
c
        dipvec_interp(1,1)=dipvec(1)
        dipvec_interp(2,1)=dipvec(2)
        dipvec_interp(3,1)=dipvec(3)

ccc        write(*,*) trinorm, dipvec
        
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         dipstr_interp(i) = dipstr*weights(i)
         call tri_bak(w,vertout,zparts(1,i))
         dipvec_interp(1,i) = dipvec_interp(1,1)
         dipvec_interp(2,i) = dipvec_interp(2,1)
         dipvec_interp(3,i) = dipvec_interp(3,1)
      enddo
      call l3dformta_dp_trunc(ier,scale,zparts,dipstr_interp,
     $   dipvec_interp,nnodes,x0y0z0,nterms,nterms,local,wlege,nlege)
      return
      end
c
c
c
c
c
      subroutine l3dformtatrid2_hess_add(ier,scale,triang,trinorm,
     $     dipstr,dipvec,
     1     ntri,x0y0z0,norder,nterms,local,wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     dipvec             = array of dipole vectors
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of local expansion
c     mptemp             = work array to hold temp local expansion
c     w(lw)              = additional work array
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale, wlege(1)
      real *8  dipvec(3,1)
      complex *16  eye, zk, dipstr(1)
      complex *16  local(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformtatridone2_hess
     $   (ier,scale,triang(1,1,i),trinorm(1,i),dipstr(i),dipvec(1,i),
     $   x0y0z0, norder,nterms,mptemp,wlege,nlege)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c
