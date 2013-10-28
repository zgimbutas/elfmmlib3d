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
c       Fortran 90 version
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a sample driver for the double layer elastostatic kernel
c     in the free space
c
c     Check accuracy for displacement and strain on a user defined grid
c
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),du(3),rnorm(3),fvec(3),strain(3,3)
        dimension triangles(3,3,1000),trinorm(3,1000),sigma(3,1000)
        dimension sigma_sl(3,1000),sigma_dl(3,1000),centroids(3,3)
        dimension fvec0(3),strain0(3,3)
c
c       SET ALL PARAMETERS
c

c     Initialize simple printing routines.
c     Calling PRINI(6,13) causes printing to screen and file fort.13.     
c
        call prini(6,13)
c
        rlam = 3.45d0
        rmu = 2.3d0
        rnu = rlam/2/(rlam+rmu)

        done=1 
        pi=4*atan(done)

        alpha = (rlam+rmu)/(rlam+2*rmu)

        x=.0d0 
        y=.0d0
        z=-1d-14

        x=-1.5d0 
        y=-1.5d0
        z=-1d-14

        x=2d0 
        y=-2d0
        z=0

        depth = 1.50d0
        dip = 0
        phi = dip/180.0d0*pi

        pot1 = -2
        pot2 = 3
        pot3 = 1
        pot4 = 0

c
        source(1)=0
        source(2)=0
        source(3)=-depth

c
c       convert strike-dip-tensile parameters into cartesian coordinates
c
        du(1)=pot1
        du(2)=pot2*cos(phi)-pot3*sin(phi)
        du(3)=pot2*sin(phi)+pot3*cos(phi)
        rnorm(1)=0
        rnorm(2)=-sin(phi)
        rnorm(3)=+cos(phi)

        scale=1d0
        du(1)=du(1)*scale
        du(2)=du(2)*scale
        du(3)=du(3)*scale

        call prin2('du=*',du,3)
        call prin2('rnorm=*',rnorm,3)

        itri=1
        triangles(1,1,itri) = 0 
        triangles(2,1,itri) = 0
        triangles(3,1,itri) = -depth
        triangles(1,2,itri) = 1
        triangles(2,2,itri) = 0
        triangles(3,2,itri) = -depth
        triangles(1,3,itri) = 1
        triangles(2,3,itri) = cos(phi)
        triangles(3,3,itri) = -depth+sin(phi)

        itri=2
        triangles(1,1,itri) = 0
        triangles(2,1,itri) = 0
        triangles(3,1,itri) = -depth
        triangles(1,2,itri) = 1
        triangles(2,2,itri) = cos(phi)
        triangles(3,2,itri) = -depth+sin(phi)
        triangles(1,3,itri) = 0
        triangles(2,3,itri) = cos(phi)
        triangles(3,3,itri) = -depth+sin(phi)

        ntri = 1

        do i = 1,ntri
        sigma(1,i)=du(1)
        sigma(2,i)=du(2)
        sigma(3,i)=du(3)
        enddo

        call triangle_norm(triangles(1,1,1),trinorm(1,1))
        call triangle_norm(triangles(1,1,2),trinorm(1,2))
        call prin2('triangles=*',triangles,3*3*ntri)
        call prin2('trinorm=*',trinorm,3*ntri)


        ifsingle=0
        ifdouble=1
        do i = 1,ntri
        sigma_sl(1,i)=0
        sigma_sl(2,i)=0
        sigma_sl(3,i)=0
        sigma_dl(1,i)=du(1)
        sigma_dl(2,i)=du(2)
        sigma_dl(3,i)=du(3)
        enddo
c
        do i = 1,ntri
        centroids(1,i)=
     $     (triangles(1,1,i)+triangles(1,2,i)+triangles(1,3,i))/3
        centroids(2,i)=
     $     (triangles(2,1,i)+triangles(2,2,i)+triangles(2,3,i))/3
        centroids(3,i)=
     $     (triangles(3,1,i)+triangles(3,2,i)+triangles(3,3,i))/3
        enddo


        write(15,*) 'triangle=['
        write(15,*) triangles(1,1,1),'',triangles(2,1,1),'',
     $     triangles(3,1,1),''
        write(15,*) triangles(1,2,1),'',triangles(2,2,1),'',
     $     triangles(3,2,1),''
        write(15,*) triangles(1,3,1),'',triangles(2,3,1),'',
     $     triangles(3,3,1),''
        write(15,*) ']'
        
        write(15,*) 'trinorm=['
        write(15,*) trinorm(1,1),'',trinorm(2,1),'',
     $     trinorm(3,1),''
        write(15,*) ']'

        write(15,*) 'sigma=['
        write(15,*) du(1),'',du(2),'',
     $     du(3),''
        write(15,*) ']'
        
        write(15,*) 'ifsingle=', ifsingle
        write(15,*) 'ifdouble=', ifdouble
        write(15,*) 'ifhalfspace=', 0


        ngridx=40
        ngridy=40
        ngridz=1
        sx=4d0
        sy=4d0
        sz=0

        ngridx=1
        ngridy=40
        ngridz=40
        sx=0d0
        sy=4d0
        sz=-3d0

        write(15,*) 'errors=['

        do 1300 ix=1,ngridx
        do 1200 iy=1,ngridy+1
        do 1100 iz=1,ngridz+1

        target(1)=x+(ix-1)/dble(ngridx)*sx
        target(2)=y+(iy-1)/dble(ngridy)*sy
        target(3)=z+(iz-1)/dble(ngridz)*sz

        ifstrain = 1


        call prin2(' *',fvec,0)
        call prin2('target=*',target,3)
        call prin2('=== Triangles (analytic) ===*',fvec,0)

        
        ntimes=1
        t1=second()
C$        t1=omp_get_wtime()
        do i=1,ntimes        
        ifptfrc=0
        ifstrain=0
        ifptfrctarg=1
        ifstraintarg=1
        ntarget=1
        call el3dtriadirecttarg(
     $     RLAM,RMU,TRIANGLES,TRINORM,NTRI,CENTROIDS,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc0,ifstrain,strain0,NTARGET,
     $     target,ifptfrctarg,fvec,
     $     ifstraintarg,STRAIN)
        enddo
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('el3dtriadirecttarg: speed, triangles/sec=*',
     $     ntimes*ntri/(t2-t1),1)

c        do i = 1,3
c        fvec(i)=fvec(i)/(4*pi)
c        enddo
c        do i = 1,3
c        do j = 1,3
c        strain(i,j)=strain(i,j)/(4*pi)
c        enddo
c        enddo

        call prin2('after el3dtriadirecttarg, fvec=*',fvec,3)
        call prin2('after el3dtriadirecttarg, strain=*',strain,3*3)

        do i=1,3
        fvec0(i)=fvec(i)
        do j=1,3
        strain0(i,j)=strain(i,j)
        enddo
        enddo

        call prin2(' *',fvec,0)
        call prin2('=== Triangles (adaptive integration) ===*',fvec,0)

        ntimes=1
        t1=second()
C$        t1=omp_get_wtime()
        do i=1,ntimes
        call elt3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,fvec,strain,numfunev)
        enddo
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('elt3triaadap: speed, triangles/sec=*',
     $     ntimes*ntri/(t2-t1),1)
        call prinf('numfunev=*',numfunev,1)

c        do i = 1,3
c        fvec(i)=fvec(i)/(4*pi)
c        enddo
c        do i = 1,3
c        do j = 1,3
c        strain(i,j)=strain(i,j)/(4*pi)
c        enddo
c        enddo

        call prin2('after elt3triaadap, fvec=*',fvec,3)
        call prin2('after elt3triaadap, strain=*',strain,3*3)
c
c       
        er1=0
        do i=1,3
        er1=er1+(fvec0(i)-fvec(i))**2
        enddo
        er1=sqrt(er1)
        call prin2('error in displacement=*',er1,1)

        er2=0
        do i=1,3
        do j=1,3
        er2=er2+(strain0(i,j)-strain(i,j))**2
        enddo
        enddo
        er2=sqrt(er2)
        call prin2('error in strain=*',er2,1)

        write(14,2000) target(1), target(2), target(3), er1, er2
        write(15,2000) target(1), target(2), target(3), er1, er2
 2000   format(6(1x,e13.5))

 1100   continue
 1200   continue
 1300   continue

        write(15,*) ']'

        stop
        end
c
c
c
c
c
        subroutine elt3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,ptfrc,strain,numfunev)
C
C     adaptive integration on triangles subroutine for elastostatic 
c     double layer N-body problem.
C
C     INPUT:
C
C     RLAM, RMU = Lame parameters
C     NTRI = number of triangles
C     TRIANGLES(3,3,ntri) = (x,y,z)-coordinate of 3 vertices 
C     SIGMA(3,ntri) = vector strength of nth charge
C     TRINORM(3,ntri) =  normal vectors
C     TARGET(3) = evaluation point
C
C     OUTPUT:
C
C     PTFRC(3) = computed displacement at target
C     STRAIN(3,3) = computed strain at target
C
        implicit real *8 (a-h,o-z)
        real *8 triangles(3,3,ntri),sigma(3,ntri)
        real *8 trinorm(3,ntri)
        real *8 target(3),source(3),ptfrc(3),strain(3,3)
        real *8 ptfrc0(3),strain0(3,3)
        real *8 vert1(3),vert2(3),vert3(3)
        dimension rints(100),par(100)
	real *8, allocatable :: w(:)
        integer ntri
c
        external fun3elt_eval
        external fun3elth_eval
c
	allocate( w(100000) )
c
        ptfrc(1) = 0.0d0
        ptfrc(2) = 0.0d0
        ptfrc(3) = 0.0d0
c
        do i=1,3
        do j=1,3
        strain(i,j) = 0.0d0
        enddo
        enddo
c       
c       ... note that the sum of weights is .5, 
c       we will need to include a factor of 2 later
c
        numfunev=0
c
        do 1200 itri = 1,ntri
c
        vert1(1) = triangles(1,1,itri)
        vert1(2) = triangles(2,1,itri)
        vert1(3) = triangles(3,1,itri)
        vert2(1) = triangles(1,2,itri)
        vert2(2) = triangles(2,2,itri)
        vert2(3) = triangles(3,2,itri)
        vert3(1) = triangles(1,3,itri)
        vert3(2) = triangles(2,3,itri)
        vert3(3) = triangles(3,3,itri)
        
        par(1)=rlam
        par(2)=rmu
        par(3)=sigma(1,itri)
        par(4)=sigma(2,itri)
        par(5)=sigma(3,itri)
        par(6)=trinorm(1,itri)
        par(7)=trinorm(2,itri)
        par(8)=trinorm(3,itri)
        
        nq = 20
        nfuns = 12
c       
        eps=1e-12
        call tria3adam(ier,vert1,vert2,vert3,fun3elt_eval,nfuns,
     1      target,par,nq,eps,rints,maxrec,numfunev0,w)
ccc        call prinf('ier=*',ier,1)
ccc        call prinf('numfunev=*',numfunev0,1)
        numfunev=numfunev+numfunev0

        k=0
        do i = 1,3
        k=k+1
        ptfrc0(i)=rints(k)
        enddo
        
        do i = 1,3
        do j = 1,3
        k=k+1
        strain0(i,j)=rints(k)
        enddo
        enddo
c
c
        do i = 1,3
        ptfrc(i) = ptfrc(i) + ptfrc0(i) 
        enddo
        
        do i = 1,3
        do j = 1,3
        strain(i,j) = strain(i,j) + strain0(i,j) 
        enddo
        enddo
c
 1200   continue
        return
        end
c
c
c
c
c
        subroutine fun3elt_eval(x,y,z,target,par,f)
        implicit real *8 (a-h,o-z)
        dimension target(3),source(3),rvec(3)
        dimension ptfrc0(3),strain0(3,3),f(1)
        dimension sigma(3),trinorm(3),par(1)
c
        source(1)=x
        source(2)=y
        source(3)=z
c
        rlam=par(1)
        rmu=par(2)
        sigma(1)=par(3)
        sigma(2)=par(4)
        sigma(3)=par(5)
        trinorm(1)=par(6)
        trinorm(2)=par(7)
        trinorm(3)=par(8)
c
        ifstrain = 1
        call green3elt_eval2
     $     (rlam,rmu,source,sigma,trinorm,target,
     $     ptfrc0,ifstrain,strain0)
c
        k=0
        do i = 1,3
        k=k+1
        f(k) = ptfrc0(i)
        enddo
        
        do i = 1,3
        do j = 1,3
        k=k+1
        f(k) = strain0(i,j)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
