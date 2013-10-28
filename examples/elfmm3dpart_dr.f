c
c       This is a simple driver to test the particle FMM routines
c       in R^3 using the free-space Green's functions for general values 
c       of the Lame parameters rlam/rmu. 
c
        implicit real *8 (a-h,o-z)
        parameter (lw=120 000 000)
        dimension w(lw)
c
        call elfmm3d_test(w,lw)
c
        stop
        end
c
c
c
c
c
        subroutine elfmm3d_test(w,lw)
        implicit real *8 (a-h,o-z)
c
c       Compute potentials by direct calculation and via FMM.
c
        parameter(nmax=100000)
c
        dimension sigma_sl(3,nmax)
        dimension sigma_dl(3,nmax)
        dimension sigma_dv(3,nmax)
c
c
        dimension source(3,nmax)
        dimension target(3,nmax)
c
        dimension w(lw),xyz(3)
c
        dimension ptfrc(3,nmax),strain(3,3,nmax)
        dimension ptfrctarg(3,nmax),straintarg(3,3,nmax)
c
        dimension ptfrc0(3),strain0(3,3)
c
        dimension ptfrc1(3,nmax),strain1(3,3,nmax)
c
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        print *, 'ENTER n'
        read *, nsource
c
c
        call prinf('nsource=*',nsource,1)
c
        idist=3
c
        if( idist .eq. 1 ) then
c
c       ... construct randomly located charge distribution on a unit cube
c       
        do i=1,nsource
        source(1,i)=hkrand(0)
        source(2,i)=hkrand(0)
        source(3,i)=hkrand(0)
        source(1,i)=source(1,i)-0.5
        source(2,i)=source(2,i)-0.5
        source(3,i)=source(3,i)-0.5
        enddo
c
        endif
c
        if( idist .eq. 2 ) then
c
c       ... construct charge distribution on a curve in R^3
c       
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        source(1,i)=sin(1.1*a)
        source(2,i)=cos(2.1*a)
        source(3,i)=cos(3.1*a)
        enddo
c
        endif
c       
        if( idist .eq. 3 ) then
c
c       ... construct randomly located charge distribution on a unit sphere
c 
        done=1
        pi=4*atan(done)
c
        d=hkrand(0)
        do i=1,nsource
c
c        source(1,i)=hkrand(0)
c        source(2,i)=hkrand(0)
c        source(3,i)=hkrand(0)
c        source(1,i)=source(1,i)-0.5
c        source(2,i)=source(2,i)-0.5
c        source(3,i)=source(3,i)-0.5
c        rr=source(1,i)**2+source(2,i)**2+source(3,i)**2
c        rr=sqrt(rr)
c        source(1,i)=source(1,i)/rr
c        source(2,i)=source(2,i)/rr
c        source(3,i)=source(3,i)/rr
c
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        source(1,i)=.5d0*cos(phi)*sin(theta)
        source(2,i)=.5d0*sin(phi)*sin(theta)
        source(3,i)=.5d0*cos(theta)
c
        d=hkrand(0)
        d=hkrand(0)
        enddo
c
        endif
c
        if( idist .eq. 4 ) then
c
c       ... construct the grid of charges
c
        ngrid=nsource-1
        kk=0
        do i=1,ngrid+1
        do j=1,ngrid+1
        do k=1,ngrid+1
        kk=kk+1
        source(1,kk)=(i-1.0)/ngrid
        source(2,kk)=(j-1.0)/ngrid
        source(3,kk)=(k-1.0)/ngrid
        source(1,kk)=source(1,kk)+hkrand(0)*.0001
        source(2,kk)=source(2,kk)+hkrand(0)*.0001
        source(3,kk)=source(3,kk)+hkrand(0)*.0001
ccc        call prin2('source=*',source(1,kk),3)
        enddo
        enddo
        enddo
        nsource=kk
c       
        call prinf('after grid build, nsource=*',nsource,1)
c
        endif
c
c        do i=1,nsource
c        source(1,i)=source(1,i)*10
c        source(2,i)=source(2,i)*10
c        source(3,i)=source(3,i)*10
c        enddo


        do i=1,nsource
        source(1,i)=source(1,i)
        source(2,i)=source(2,i)
        source(3,i)=source(3,i) - 3
        enddo

c       ... set up the targets
c
        if( idist .eq. 1 .or. idist .eq. 4 ) then
        do i=1,nsource
        target(1,i)=source(1,i) + 0.1
        target(2,i)=source(2,i)
        target(3,i)=source(3,i)
        enddo
        ntarget=nsource
        do i=1,nsource
        target(1,i+nsource)=source(1,i) 
        target(2,i+nsource)=source(2,i) - 0.2
        target(3,i+nsource)=source(3,i)
        enddo
        ntarget=nsource*2
        endif
c
c
        if( idist .eq. 2 ) then
c
c       ... construct target distribution on a curve in R^3
c       
        ntarget=nsource*4
        do i=1,ntarget
        a=2*pi*dble(i)/ntarget
        target(1,i)=sin(1.1*a)/2
        target(2,i)=cos(2.1*a)/2
        target(3,i)=cos(3.1*a)/2
        enddo
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct target distribution on a unit sphere 
c       highly oversampled 
c
        ntarget=nsource*1
        do i=1,ntarget
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        target(1,i)=.5d0*cos(phi)*sin(theta) + 1
        target(2,i)=.5d0*sin(phi)*sin(theta)
        target(3,i)=.5d0*cos(theta) - 2
        d=hkrand(0)
        d=hkrand(0)
        enddo
        endif
c
        call prinf('ntarget=*',ntarget,1)
c       
ccc        call prin2('source=*',source,3*nsource)
ccc        call prin2('target=*',target,3*ntarget)
        
c
c       define (random) piecewise constant densities
c
        do i=1,nsource
        sigma_sl(1,i)=hkrand(0)
        sigma_sl(2,i)=hkrand(0)
        sigma_sl(3,i)=hkrand(0)
        sigma_dl(1,i)=hkrand(0)
        sigma_dl(2,i)=hkrand(0)
        sigma_dl(3,i)=hkrand(0)
        enddo        
c
        do i=1,nsource
        sigma_dv(1,i)=hkrand(0)
        sigma_dv(2,i)=hkrand(0)
        sigma_dv(3,i)=hkrand(0)
        enddo        
c
c
c       
        rlam=3.2d0
        rmu=2.3d0
c
c       turn on single and/or double layer potential
c
        ifsingle=1
        ifdouble=1
c
c       set whether displacement and/or strain to be compute on surface
c       and at target locations.
c
        ifptfrc=1
        ifstrain=1
        ifptfrctarg=1
        ifstraintarg=1
c
        ifprint=0
c
c       ... evaluate via FMM 
c
        iprec=1
c
ccc        call prini(0,13)
        t1=second()
C$        t1=omp_get_wtime()
        call elfmm3dparttarg
     $     (ier,iprec,RLAM,RMU,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prini(6,13)
c
        call prin2('after elfmm3parttarg, time=*',t2-t1,1)
        call prin2('speed, (sources+targets)/sec=*',
     $     (nsource+ntarget)/(t2-t1),1)
c
c



c       since we are testing by direct calculation, only compute
c       direct calculation at m sources.
c
        m=min(nsource,10)
c
c       ... direct evaluation
c
        do i=1,m
        do j=1,3
           ptfrc1(j,i)=0
        enddo
        do j=1,3
        do k=1,3
           strain1(j,k,i)=0
        enddo
        enddo
        enddo
c
        if( ifprint .eq. 1 ) then
        call prin2('after elfmm3parttarg, ptfrc=*',
     $     ptfrc,3*m)
        call prin2('after elfmm3parttarg, strain=*',
     $     strain,3*3*m)
        endif
c
c       
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do 6550 j=1,m
        do 6540 i=1,nsource
        if( i .eq. j ) goto 6540

        if (ifsingle .eq. 1 ) then
        call green3elu_eval(rlam,rmu,source(1,i),
     $     sigma_sl(1,i),source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc1(1,j)=ptfrc1(1,j)+ptfrc0(1)
        ptfrc1(2,j)=ptfrc1(2,j)+ptfrc0(2)
        ptfrc1(3,j)=ptfrc1(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain1(1,1,j)=strain1(1,1,j)+strain0(1,1)
        strain1(2,1,j)=strain1(2,1,j)+strain0(2,1)
        strain1(3,1,j)=strain1(3,1,j)+strain0(3,1)
        strain1(1,2,j)=strain1(1,2,j)+strain0(1,2)
        strain1(2,2,j)=strain1(2,2,j)+strain0(2,2)
        strain1(3,2,j)=strain1(3,2,j)+strain0(3,2)
        strain1(1,3,j)=strain1(1,3,j)+strain0(1,3)
        strain1(2,3,j)=strain1(2,3,j)+strain0(2,3)
        strain1(3,3,j)=strain1(3,3,j)+strain0(3,3)
        endif
        endif
c
        if (ifdouble .eq. 1) then
        call green3elt_eval(rlam,rmu,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),ptfrc0,ifstrain,strain0)
        if (ifptfrc .eq. 1) then
        ptfrc1(1,j)=ptfrc1(1,j)+ptfrc0(1)
        ptfrc1(2,j)=ptfrc1(2,j)+ptfrc0(2)
        ptfrc1(3,j)=ptfrc1(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain1(1,1,j)=strain1(1,1,j)+strain0(1,1)
        strain1(2,1,j)=strain1(2,1,j)+strain0(2,1)
        strain1(3,1,j)=strain1(3,1,j)+strain0(3,1)
        strain1(1,2,j)=strain1(1,2,j)+strain0(1,2)
        strain1(2,2,j)=strain1(2,2,j)+strain0(2,2)
        strain1(3,2,j)=strain1(3,2,j)+strain0(3,2)
        strain1(1,3,j)=strain1(1,3,j)+strain0(1,3)
        strain1(2,3,j)=strain1(2,3,j)+strain0(2,3)
        strain1(3,3,j)=strain1(3,3,j)+strain0(3,3)
        endif
        endif
 6540   continue
 6550   continue
C$OMP END PARALLEL DO

        if( ifprint .eq. 1 ) then
        call prin2('after directtarg, ptfrc=*',ptfrc1,3*m)
        call prin2('after directtarg, strain=*',strain1,3*3*m)
        endif
c
        t2=second()
C$        t2=omp_get_wtime()

        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)

        if( ifptfrc .eq. 1 ) then
        call d3error(ptfrc1,ptfrc,3*m,a,r)
ccc        call prin2('absolute error in ptfrc=*',a,1)
        call prin2('relative error in ptfrc=*',r,1)
        endif
c
        if( ifstrain .eq. 1 ) then
        call d3error(strain1,strain,3*3*m,a,r)
ccc        call prin2('absolute error in strain=*',a,1)
        call prin2('relative error in strain=*',r,1)
        endif
c
c

        if( ifptfrctarg .eq. 0 .and. ifstraintarg .eq. 0 ) return


c       since we are testing by direct calculation, only compute
c       direct calculation at m targets.
c
        m=min(ntarget,10)
c

c       ... direct evaluation
c
        do i=1,m
        do j=1,3
           ptfrc1(j,i)=0
        enddo
        do j=1,3
        do k=1,3
           strain1(j,k,i)=0
        enddo
        enddo
        enddo
c
        if( ifprint .eq. 1 ) then
        call prin2('after elfmm3parttarg, ptfrctarg=*',
     $     ptfrctarg,3*m)
        call prin2('after elfmm3parttarg, straintarg=*',
     $     straintarg,3*3*m)
        endif
c

c       
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0)
        do j=1,m
        do i=1,nsource
        if (ifsingle .eq. 1 ) then
        call green3elu_eval(rlam,rmu,source(1,i),
     $     sigma_sl(1,i),target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrc1(1,j)=ptfrc1(1,j)+ptfrc0(1)
        ptfrc1(2,j)=ptfrc1(2,j)+ptfrc0(2)
        ptfrc1(3,j)=ptfrc1(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        strain1(1,1,j)=strain1(1,1,j)+strain0(1,1)
        strain1(2,1,j)=strain1(2,1,j)+strain0(2,1)
        strain1(3,1,j)=strain1(3,1,j)+strain0(3,1)
        strain1(1,2,j)=strain1(1,2,j)+strain0(1,2)
        strain1(2,2,j)=strain1(2,2,j)+strain0(2,2)
        strain1(3,2,j)=strain1(3,2,j)+strain0(3,2)
        strain1(1,3,j)=strain1(1,3,j)+strain0(1,3)
        strain1(2,3,j)=strain1(2,3,j)+strain0(2,3)
        strain1(3,3,j)=strain1(3,3,j)+strain0(3,3)
        endif
        endif
c
        if (ifdouble .eq. 1) then
        call green3elt_eval(rlam,rmu,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrc1(1,j)=ptfrc1(1,j)+ptfrc0(1)
        ptfrc1(2,j)=ptfrc1(2,j)+ptfrc0(2)
        ptfrc1(3,j)=ptfrc1(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        strain1(1,1,j)=strain1(1,1,j)+strain0(1,1)
        strain1(2,1,j)=strain1(2,1,j)+strain0(2,1)
        strain1(3,1,j)=strain1(3,1,j)+strain0(3,1)
        strain1(1,2,j)=strain1(1,2,j)+strain0(1,2)
        strain1(2,2,j)=strain1(2,2,j)+strain0(2,2)
        strain1(3,2,j)=strain1(3,2,j)+strain0(3,2)
        strain1(1,3,j)=strain1(1,3,j)+strain0(1,3)
        strain1(2,3,j)=strain1(2,3,j)+strain0(2,3)
        strain1(3,3,j)=strain1(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
        enddo
C$OMP END PARALLEL DO

        t2=second()
C$        t2=omp_get_wtime()

        if( ifprint .eq. 1 ) then
        call prin2('after directtarg, ptfrctarg=*',ptfrc1,3*m)
        call prin2('after directtarg, straintarg=*',strain1,3*3*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)

        if( ifptfrctarg .eq. 1 ) then
        call d3error(ptfrc1,ptfrctarg,3*m,a,r)
ccc        call prin2('absolute error in target ptfrc=*',a,1)
        call prin2('relative error in target ptfrc=*',r,1)
        endif
c
        if( ifstraintarg .eq. 1 ) then
        call d3error(strain1,straintarg,3*3*m,a,r)
ccc        call prin2('absolute error in target strain=*',a,1)
        call prin2('relative error in target strain=*',r,1)
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine d3error(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        dimension pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
        d=d+abs(pot1(i)-pot2(i))**2
        a=a+abs(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
