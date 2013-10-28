cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-01-05 16:30:53 -0500 (Tue, 05 Jan 2010) $
c    $Revision: 786 $
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains the FMM routines for elastostatic particle
c        potentials in free space in R^3. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      elfmm3dpart - Elastostatic FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c      elfmm3dparttarg - Elastostatic FMM in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c      el3dpartdirecttarg - Elastostatic interactions in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
c
c      elh3dpartdirecttarg - Elastostatic interactions in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm. Half space
c         elastostatic Green's function (Mindlin's solution).
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine elfmm3dpart
     $     (ier,iprec,rlam,rmu,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain)
c
c     FMM calculation subroutine for elastostatic N-body problem.
c
c     Free space elastostatic Green's function.
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     ptfrc(3,nparts) = displacement at source locations
c     strain(3,3,nparts) = strain at source locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        integer nparts,ntargs

        ntargs=0
        ifptfrctarg=0
        ifstraintarg=0
        call elfmm3dparttarg
     $     (ier,iprec,RLAM,RMU,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)

        return
        end
c
c
c
c
c*********************************
      subroutine elfmm3dparttarg
     $     (ier,iprec,rlam,rmu,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
c
c     FMM calculation subroutine for elastostatic N-body problem.
c
c     Free space elastostatic Green's function.
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     ptfrc(3,nparts) = displacement 
c     strain(3,3,nparts) = strain
c     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrctarg(3,ntargs),straintarg(3,3,ntargs)
        integer nparts,ntargs
        real *8, allocatable :: w(:)
c       
        ier=0
        lused=0
c
c       ... allocate work arrays
c
        icharge=lused+1
        lcharge=2*nparts
        lused=lused+lcharge

        idipstr=lused+1
        ldipstr=2*nparts
        lused=lused+ldipstr

        idipvec=lused+1
        ldipvec=3*nparts
        lused=lused+ldipvec

        ipot=lused+1
        lpot=2*nparts
        lused=lused+lpot

        ifld=lused+1
        lfld=2*3*nparts
        lused=lused+lfld

        ihess=lused+1
        lhess=2*6*nparts
        lused=lused+lhess

        ihessmatr=lused+1
        lhessmatr=3*3*nparts
        lused=lused+lhessmatr

        ipottarg=lused+1
        lpottarg=2*ntargs
        lused=lused+lpottarg

        ifldtarg=lused+1
        lfldtarg=2*3*ntargs
        lused=lused+lfldtarg

        ihesstarg=lused+1
        lhesstarg=2*6*ntargs
        lused=lused+lhesstarg

        ihessmatrtarg=lused+1
        lhessmatrtarg=3*3*ntargs
        lused=lused+lhessmatrtarg
c
        allocate(w(lused+10),stat=ier)
        if( ier .ne. 0 ) return
c
        call elfmm3dparttargmain_fast
     $     (ier,iprec,RLAM,RMU,NPARTS,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg,
     $     w(icharge),w(idipstr),w(idipvec),
     $     w(ipot),w(ifld),w(ihess),w(ihessmatr),
     $     w(ipottarg),w(ifldtarg),w(ihesstarg),w(ihessmatrtarg))
c
        if( ier .ne. 0 ) return
c
        return
        end
c
c
c
c
c
C*********************************
      subroutine elfmm3dparttargmain_fast
     $     (ier,iprec,rlam,rmu,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg,
     $     charge,dipstr,dipvec,pot,fld,hess,hessmatr,
     $     pottarg,fldtarg,hesstarg,hessmatrtarg)
c
c     FMM calculation subroutine for elastostatic N-body problem
c
c     4 Laplace FMM calls.
c
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     ptfrc(3,nparts) = displacement
c     strain(3,3,nparts) = strain
c     ptfrctarg(3,ntargs) = displacement at target locations
c     straintarg(3,3,ntargs) = strain at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 ptfrc(3,nparts),strain(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 ptfrctarg(3,ntargs),straintarg(3,3,ntargs)
        integer nparts,ntargs
c       
        complex *16 charge(1)
        complex *16 dipstr(1)
        real *8 dipvec(3,1)
c
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1) 
        real *8 hessmatr(3,3,1)

        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        complex *16 hesstarg(6,1) 
        real *8 hessmatrtarg(3,3,1)
c
c
        do k=1,nparts
c
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,k) = 0.0d0
        ptfrc(2,k) = 0.0d0
        ptfrc(3,k) = 0.0d0
        endif
c
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        hessmatr(i,j,k) = 0.0d0
        enddo
        enddo
        endif        
c
        enddo
c
c
        do k=1,ntargs
c
        if( ifptfrctarg .eq. 1 ) then
        ptfrctarg(1,k) = 0.0d0
        ptfrctarg(2,k) = 0.0d0
        ptfrctarg(3,k) = 0.0d0
        endif
c
        if( ifstraintarg .eq. 1 ) then
        do i=1,3
        do j=1,3
        hessmatrtarg(i,j,k) = 0.0d0
        enddo
        enddo
        endif
c
        enddo
c
c
c
        ifpot=0
        iffld=0
        ifhess=0
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
        ifpot=1
        iffld=1
        endif
        if( ifstrain .eq. 1 ) ifhess=1

        ifpottarg=0
        iffldtarg=0
        ifhesstarg=0
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
        ifpottarg=1
        iffldtarg=1
        endif
        if( ifstraintarg .eq. 1 ) ifhesstarg=1
c
c
        c1 = (rlam+rmu)/(rlam+2*rmu)
        c2 = (rmu)/(rlam+2*rmu)

c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
        do j = 1,3

        ifcharge=0
        ifdipole=0

        do k = 1,nparts
            charge(k) = 0
            dipstr(k) = 0
            dipvec(1,k) = 0
            dipvec(2,k) = 0
            dipvec(3,k) = 0
            if( ifsingle .eq. 1 ) then
            charge(k) = sigma_sl(j,k)/(2*rmu)
            ifcharge=1
            endif
            if( ifdouble .eq. 1 ) then
            dipstr(k) = 1
            dipvec(1,k) = sigma_dv(1,k)*sigma_dl(j,k)
            dipvec(2,k) = sigma_dv(2,k)*sigma_dl(j,k)
            dipvec(3,k) = sigma_dv(3,k)*sigma_dl(j,k)
            dipvec(1,k) = dipvec(1,k)+sigma_dl(1,k)*sigma_dv(j,k)
            dipvec(2,k) = dipvec(2,k)+sigma_dl(2,k)*sigma_dv(j,k)
            dipvec(3,k) = dipvec(3,k)+sigma_dl(3,k)*sigma_dv(j,k)
            dipvec(1,k) = dipvec(1,k)/2
            dipvec(2,k) = dipvec(2,k)/2
            dipvec(3,k) = dipvec(3,k)/2
            ifdipole=1
            endif
        enddo

        call lfmm3dparthesstarg(ier,iprec,
     $     nparts,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntargs,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)

        call elfmm3dlap1(nparts,j,c1,c2,pot,fld,hess,
     $     source,ifptfrc,ptfrc,ifstrain,hessmatr)
        call elfmm3dlap1(ntargs,j,c1,c2,pottarg,fldtarg,hesstarg,
     $     target,ifptfrctarg,ptfrctarg,ifstraintarg,hessmatrtarg)

        enddo

c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
        ifcharge=0
        ifdipole=0
c        
        do k = 1,nparts
          charge(k) = 0
          dipstr(k) = 0
          dipvec(1,k) = 0
          dipvec(2,k) = 0
          dipvec(3,k) = 0
          if( ifsingle .eq. 1 ) then
          charge(k) = 
     $      (sigma_sl(1,k)*source(1,k)+
     $       sigma_sl(2,k)*source(2,k)+
     $       sigma_sl(3,k)*source(3,k))*c1/(2*rmu)
          ifcharge = 1
          endif
          if( ifdouble .eq. 1 ) then
          charge(k) = charge(k) + 
     $        (sigma_dl(1,k)*sigma_dv(1,k) + 
     1         sigma_dl(2,k)*sigma_dv(2,k) + 
     2         sigma_dl(3,k)*sigma_dv(3,k))*c2
          dipstr(k) = 1
          dipvec(1,k) = c1*sigma_dv(1,k)*
     $        (sigma_dl(1,k)*source(1,k) + 
     1         sigma_dl(2,k)*source(2,k) + 
     2         sigma_dl(3,k)*source(3,k) )
          dipvec(2,k) = c1*sigma_dv(2,k)*
     $        (sigma_dl(1,k)*source(1,k) + 
     1         sigma_dl(2,k)*source(2,k) + 
     2         sigma_dl(3,k)*source(3,k) )
          dipvec(3,k) = c1*sigma_dv(3,k)*
     $        (sigma_dl(1,k)*source(1,k) + 
     1         sigma_dl(2,k)*source(2,k) + 
     2         sigma_dl(3,k)*source(3,k) )
          dipvec(1,k) = dipvec(1,k) + c1*sigma_dl(1,k)*
     $        (sigma_dv(1,k)*source(1,k) + 
     1         sigma_dv(2,k)*source(2,k) + 
     2         sigma_dv(3,k)*source(3,k))
          dipvec(2,k) = dipvec(2,k) + c1*sigma_dl(2,k)*
     $        (sigma_dv(1,k)*source(1,k) + 
     1         sigma_dv(2,k)*source(2,k) + 
     2         sigma_dv(3,k)*source(3,k))
          dipvec(3,k) = dipvec(3,k) + c1*sigma_dl(3,k)*
     $        (sigma_dv(1,k)*source(1,k) + 
     1         sigma_dv(2,k)*source(2,k) + 
     2         sigma_dv(3,k)*source(3,k))
          dipvec(1,k) = dipvec(1,k)/2
          dipvec(2,k) = dipvec(2,k)/2
          dipvec(3,k) = dipvec(3,k)/2
          ifcharge = 1
          ifdipole = 1
          endif
        enddo

        call lfmm3dparthesstarg(ier,iprec,
     $     nparts,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntargs,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)

        call elfmm3dlap2(nparts,pot,fld,hess,
     $     ifptfrc,ptfrc,ifstrain,hessmatr)
        call elfmm3dlap2(ntargs,pottarg,fldtarg,hesstarg,
     $     ifptfrctarg,ptfrctarg,ifstraintarg,hessmatrtarg)


        do k=1,nparts
c       
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,k)=ptfrc(1,k)
        ptfrc(2,k)=ptfrc(2,k)
        ptfrc(3,k)=ptfrc(3,k)
        endif
c       
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j,k)=(hessmatr(i,j,k)+hessmatr(j,i,k))/2
        enddo
        enddo
        endif
c       
        enddo
c 
        do k=1,ntargs
c       
        if( ifptfrctarg .eq. 1 ) then
        ptfrctarg(1,k)=ptfrctarg(1,k)
        ptfrctarg(2,k)=ptfrctarg(2,k)
        ptfrctarg(3,k)=ptfrctarg(3,k)
        endif
c       
        if( ifstraintarg .eq. 1 ) then
        do i=1,3
        do j=1,3
        straintarg(i,j,k)=(hessmatrtarg(i,j,k)+hessmatrtarg(j,i,k))/2
        enddo
        enddo
        endif
c
        enddo
c       
        return
        end
c
c
c
c
c
        subroutine el3dpartdirecttarg(
     $     rlam,rmu,nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,ntarget,
     $     target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
        implicit real *8 (a-h,o-z)
c
c
c       Elastostatic interactions in R^3: evaluate all pairwise particle
c       interactions (excluding self interactions) and interactions with
c       targets using the direct O(N^2) algorithm.
c
c       INPUT:
c
c       rlam,rmu - Lame parameters
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       sigma_dv(3,nsource) - orientation of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifptfrc - displacement computation flag
c       ifstrain - strain computation flag
c       ifptfrctarg - target displacement computation flag
c       ifstraintarg - target strain computation flag
c
c       OUTPUT:
c
c       ptfrc(3,nsource) - displacement at source locations
c       strain(3,3,nsource) - strain at source locations
c       ptfrctarg(3,ntarget) - displacement at target locations
c       straintarg(3,3,ntarget) - strain at target locations
c
        dimension source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
c
c
        do i=1,nsource
        if( ifptfrc .eq. 1) then
           ptfrc(1,i)=0
           ptfrc(2,i)=0
           ptfrc(3,i)=0
        endif
        if( ifstrain .eq. 1) then
           strain(1,1,i)=0
           strain(2,1,i)=0
           strain(3,1,i)=0
           strain(1,2,i)=0
           strain(2,2,i)=0
           strain(3,2,i)=0
           strain(1,3,i)=0
           strain(2,3,i)=0
           strain(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifptfrctarg .eq. 1) then
           ptfrctarg(1,i)=0
           ptfrctarg(2,i)=0
           ptfrctarg(3,i)=0
        endif
        if( ifstraintarg .eq. 1) then
           straintarg(1,1,i)=0
           straintarg(2,1,i)=0
           straintarg(3,1,i)=0
           straintarg(1,2,i)=0
           straintarg(2,2,i)=0
           straintarg(3,2,i)=0
           straintarg(1,3,i)=0
           straintarg(2,3,i)=0
           straintarg(3,3,i)=0
        endif
        enddo
c
c       ... sources
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do 6550 j=1,nsource
        do 6540 i=1,nsource
        if( i .eq. j ) goto 6540

        if (ifsingle .eq. 1 ) then
        call green3elu_eval(rlam,rmu,source(1,i),
     $     sigma_sl(1,i),source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif

        if (ifdouble .eq. 1) then
        call green3elt_eval(rlam,rmu,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
c
 6540   continue
 6550   continue
c
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do j=1,ntarget
        do i=1,nsource
        if (ifsingle .eq. 1 ) then
        call green3elu_eval(rlam,rmu,source(1,i),
     $     sigma_sl(1,i),target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
c
        if (ifdouble .eq. 1) then
        call green3elt_eval(rlam,rmu,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
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
c
        subroutine el3dpartdirecttargtime(
     $     rlam,rmu,ms,nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifptfrc,ptfrc,ifstrain,strain,mt,ntarget,
     $     target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)
        implicit real *8 (a-h,o-z)
c
c
c       Elastostatic interactions in R^3: evaluate all pairwise particle
c       interactions (excluding self interactions) and interactions with
c       targets using the direct O(N^2) algorithm.
c
c       INPUT:
c
c       rlam,rmu - Lame parameters
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       sigma_dv(3,nsource) - orientation of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifptfrc - displacement computation flag
c       ifstrain - strain computation flag
c       ifptfrctarg - target displacement computation flag
c       ifstraintarg - target strain computation flag
c
c       OUTPUT:
c
c       ptfrc(3,nsource) - displacement at source locations
c       strain(3,3,nsource) - strain at source locations
c       ptfrctarg(3,ntarget) - displacement at target locations
c       straintarg(3,3,ntarget) - strain at target locations
c
        dimension source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
c
c
        do i=1,nsource
        if( ifptfrc .eq. 1) then
           ptfrc(1,i)=0
           ptfrc(2,i)=0
           ptfrc(3,i)=0
        endif
        if( ifstrain .eq. 1) then
           strain(1,1,i)=0
           strain(2,1,i)=0
           strain(3,1,i)=0
           strain(1,2,i)=0
           strain(2,2,i)=0
           strain(3,2,i)=0
           strain(1,3,i)=0
           strain(2,3,i)=0
           strain(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifptfrctarg .eq. 1) then
           ptfrctarg(1,i)=0
           ptfrctarg(2,i)=0
           ptfrctarg(3,i)=0
        endif
        if( ifstraintarg .eq. 1) then
           straintarg(1,1,i)=0
           straintarg(2,1,i)=0
           straintarg(3,1,i)=0
           straintarg(1,2,i)=0
           straintarg(2,2,i)=0
           straintarg(3,2,i)=0
           straintarg(1,3,i)=0
           straintarg(2,3,i)=0
           straintarg(3,3,i)=0
        endif
        enddo
c
c       ... sources
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do 6550 j=1,ms
        do 6540 i=1,nsource
        if( i .eq. j ) goto 6540

        if (ifsingle .eq. 1 ) then
        call green3elu_eval(rlam,rmu,source(1,i),
     $     sigma_sl(1,i),source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif

        if (ifdouble .eq. 1) then
        call green3elt_eval(rlam,rmu,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
c
 6540   continue
 6550   continue
c
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do j=1,mt
        do i=1,nsource
        if (ifsingle .eq. 1 ) then
        call green3elu_eval(rlam,rmu,source(1,i),
     $     sigma_sl(1,i),target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
c
        if (ifdouble .eq. 1) then
        call green3elt_eval(rlam,rmu,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
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
c
