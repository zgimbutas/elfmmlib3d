c
c       This is a simple driver to test the layer potential FMM routines
c       in R^3 using the free-space Green's functions for general values 
c       of the Lame parameters rlam/rmu. The geometry is assumed to
c       consist of a collection of flat triangles with piecewise constant
c       source densities (tractions for the single layer and/or 
c       jumps in displacement for the double layer).
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
c
c       Compute layer potentials by direct calculation and via FMM.
c
c
        implicit real *8 (a-h,o-z)
        parameter(nmax=100000)
c
        real *8 triangles(3,3,nmax),centroids(3,nmax)
        real *8 trinorm(3,nmax),triarea(nmax)
        real *8 verts(3,nmax),ifaces(3,nmax)       
        dimension itrivert(3,nmax)
c
        real *8 sigma_sl(3,nmax)
        real *8 sigma_dl(3,nmax)
c
        real *8 target(3,nmax)
c
        dimension vert1(3),vert2(3),vert3(3),vertout(3)
        dimension w(lw)
c
        dimension ptfrc(3,nmax),strain(3,3,nmax)
        dimension ptfrctarg(3,nmax),straintarg(3,3,nmax)
c
        dimension ptfrc0(3),strain0(3,3),stress0(3,3),tract0(3)
c
        dimension ptfrc1(3,nmax),strain1(3,3,nmax)
c
c
c       SET ALL PARAMETERS
c
c
c       PRINTING: prini determines whether output is printed or not
c       with subsequent calls to prin2 or prinf. prina(i1,i2) simply
c       prints to the Fortran unit numbers i1 and i2 assuming they are
c       nonzero. Setting i1 or i2 to zero suppresses printing. 
c
        call prini(6,13)
c
c       ... get scatterer geometry
c
        call getgeom(ntri,triangles,centroids,trinorm,triarea,nmax,
     1     verts,nmax,itrivert,nmax,iergeom)
        call prinf('geometry error flag is *',iergeom,1)
c
        call prinf('ntri=*',ntri,1)
c
ccc        call prin2('triangles=*',triangles,3*3*ntri)
ccc        call prin2('centroids=*',centroids,3*ntri)
ccc        call prin2('trinorm=*',trinorm,3*ntri)
c       
c
c       define (random) piecewise constant densities
c
        do i=1,ntri
        sigma_sl(1,i)=hkrand(0)
        sigma_sl(2,i)=hkrand(0)
        sigma_sl(3,i)=hkrand(0)
        sigma_dl(1,i)=hkrand(0)
        sigma_dl(2,i)=hkrand(0)
        sigma_dl(3,i)=hkrand(0)
        enddo        
c
c       define piecewise constant densities
c
        do i=1,ntri
c        sigma_sl(1,i)=1
c        sigma_sl(2,i)=1
c        sigma_sl(3,i)=1
c        sigma_dl(1,i)=1
c        sigma_dl(2,i)=1
c        sigma_dl(3,i)=1
        enddo        
c
c       set targets on one surface just outside the sphere and
c       one surface just inside the sphere.
c
        do i=1,ntri
        h=1d-4
        target(1,2*i-1)=centroids(1,i)+h*trinorm(1,i) 
        target(2,2*i-1)=centroids(2,i)+h*trinorm(2,i) 
        target(3,2*i-1)=centroids(3,i)+h*trinorm(3,i) 
        target(1,2*i  )=centroids(1,i)-h*trinorm(1,i)
        target(2,2*i  )=centroids(2,i)-h*trinorm(2,i)
        target(3,2*i  )=centroids(3,i)-h*trinorm(3,i)
        enddo
        ntargs=2*ntri
c        ntargs=0
c
        call prinf('ntargs=*',ntargs,1)
ccc        call prin2('target=*',target,3*ntargs)
c
        call prinf('==== targets ====*',i,0)
c       
c       since we are testing by direct calculation, only compute
c       direct calculation at m targets.
c
        m=min(ntargs,20)
        m=20
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
c
        call prinf('m=*',m,1)
c
c       ... evaluate via direct quadrature on triangles
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
        t1=second()
C$        t1=omp_get_wtime()
        if( ifsingle .eq. 1 ) then
           do i=1,m
              call elust3triadirecttarg
     $           (rlam,rmu,ntri,triangles,sigma_sl,
     1           target(1,i),ptfrc0,ifstraintarg,strain0)
              ptfrc1(1,i)=ptfrc1(1,i)+ptfrc0(1)
              ptfrc1(2,i)=ptfrc1(2,i)+ptfrc0(2)
              ptfrc1(3,i)=ptfrc1(3,i)+ptfrc0(3)
              if( ifstraintarg .eq. 1 ) then
              strain1(1,1,i)=strain1(1,1,i)+strain0(1,1)
              strain1(2,1,i)=strain1(2,1,i)+strain0(2,1)
              strain1(3,1,i)=strain1(3,1,i)+strain0(3,1)
              strain1(1,2,i)=strain1(1,2,i)+strain0(1,2)
              strain1(2,2,i)=strain1(2,2,i)+strain0(2,2)
              strain1(3,2,i)=strain1(3,2,i)+strain0(3,2)
              strain1(1,3,i)=strain1(1,3,i)+strain0(1,3)
              strain1(2,3,i)=strain1(2,3,i)+strain0(2,3)
              strain1(3,3,i)=strain1(3,3,i)+strain0(3,3)
              endif
           enddo
        endif
c
        if( ifdouble .eq. 1 ) then
           do i=1,m
              call eltst3triadirecttarg
     $           (rlam,rmu,ntri,triangles,sigma_dl,trinorm,
     $           target(1,i),ptfrc0,ifstraintarg,strain0)
              ptfrc1(1,i)=ptfrc1(1,i)+ptfrc0(1)
              ptfrc1(2,i)=ptfrc1(2,i)+ptfrc0(2)
              ptfrc1(3,i)=ptfrc1(3,i)+ptfrc0(3)
              if( ifstraintarg .eq. 1 ) then
              strain1(1,1,i)=strain1(1,1,i)+strain0(1,1)
              strain1(2,1,i)=strain1(2,1,i)+strain0(2,1)
              strain1(3,1,i)=strain1(3,1,i)+strain0(3,1)
              strain1(1,2,i)=strain1(1,2,i)+strain0(1,2)
              strain1(2,2,i)=strain1(2,2,i)+strain0(2,2)
              strain1(3,2,i)=strain1(3,2,i)+strain0(3,2)
              strain1(1,3,i)=strain1(1,3,i)+strain0(1,3)
              strain1(2,3,i)=strain1(2,3,i)+strain0(2,3)
              strain1(3,3,i)=strain1(3,3,i)+strain0(3,3)
              endif
           enddo
        endif
c
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after elfmm3triadirect, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (m)/(t2-t1),1)
        call prin2('after estimated time for direct=*',
     $     ntri/dble(m)*(t2-t1),1)
c
c       ... evaluate via FMM on triangles 
c
        iprec=1
c
ccc        call prini(0,13)
        t1=second()
C$        t1=omp_get_wtime()
        call elfmm3dtriatarg
     $     (ier,iprec,RLAM,RMU,TRIANGLES,TRINORM,NTRI,centroids,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntargs,target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prini(6,13)
c
        call prin2('after elfmm3triatarg, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (ntargs)/(t2-t1),1)
c
c
        call prin2('after elfmm3dtriatarg, ptfrc=*',ptfrc,3*m)
        call prin2('after elfmm3dtriatarg, strain=*',strain,3*3*m)
c
        call prin2('after triadirecttarg, ptfrc=*',ptfrc1,3*m)
        call prin2('after triadirecttarg, strain=*',strain1,3*3*m)
c
        if( ifptfrctarg .eq. 1 )
     $     call prin2('after elfmm3triatarg, ptfrctarg=*',
     $     ptfrctarg,3*m)
        if( ifstraintarg .eq. 1 )
     $     call prin2('after elfmm3triatarg, straintarg=*',
     $     straintarg,3*3*m)
c

        call d3error(ptfrc1,ptfrctarg,3*m,a,r)
        call prin2('absolute error in ptfrc=*',a,1)
        call prin2('relative error in ptfrc=*',r,1)
c
        call d3error(strain1,straintarg,3*3*m,a,r)
        call prin2('absolute error in strain=*',a,1)
        call prin2('relative error in strain=*',r,1)
c
c
        call prin2('after elfmm3triatarg, self ptfrc=*',
     $     ptfrc,3)
        call prin2('after elfmm3triatarg, self strain=*',
     $     strain,3*3)
c

        return
        end
c
c
c
c
c
        subroutine getgeom(ntri,triangles,centroids,trinorm,triarea,
     1     nmax,verts,lv,itrivert,li,iergeom)
        implicit real *8 (a-h,o-z)
        dimension triangles(3,3,1),centroids(3,1)
        dimension trinorm(3,1),triarea(1)
        dimension verts(3,lv),itrivert(3,li)
c       
c       INPUT 
c       
c       nmax = dimension of trisngles, centroids, trinorm, triarea arrays
c       verts =  work array to read in vertices 
c       lv = dimension of verts (work) array
c       itrivert = integer work array to read in vertices associated 
c          with each triangle
c       li = dimension of itrivert (work) array
c       
c       OUTPUT 
c
c       ntri = number of triangles
c       triangles = triangles array
c       centroids = centroids array
c       trinorm = triangle normals array
c       triarea = triangle areas array
c       iergeom = error return code
c      
c       iergeom = 0  normal execution 
c       iergeom = 1  error opening unit ir in atrireadchk
c       iergeom = 2  lv is too small to read all vertices
c       iergeom = 3  li is too small (must exceed ntri)
c       iergeom = 4  nmax is too small (must exceed ntri)
c
c-----------------------------------------------------------------------
c
        ir = 17
        open (unit = ir,file='sphere2880.a.tri')
ccc        open (unit = ir,file='../../Data/sphere2880.a.tri')
ccc        open (unit = ir,file='sphere20480.a.tri')
ccc        open (unit = ir,file='sphere180.a.tri')
cc        open (unit = ir,file='../../Data/sphere11520.a.tri')
cc        open (unit = ir,file='../../Data/sphere20480.a.tri')
        call atrireadchk(ir,verts,lv,nverts,itrivert,li,ntri,iergeom)
c
        if (ntri.gt.nmax) then
        iergeom = 4
        return
        endif
c    
c       create triangles from tri format data
c
c
c       scale and offset if desired.
c
        scale=1
        dx=0
        dy=0
        dz=0
c       
        do itri = 1,ntri
        v11 = verts(1,itrivert(1,itri)) *scale + dx
        v21 = verts(2,itrivert(1,itri)) *scale + dy
        v31 = verts(3,itrivert(1,itri)) *scale + dz
        v12 = verts(1,itrivert(2,itri)) *scale + dx 
        v22 = verts(2,itrivert(2,itri)) *scale + dy
        v32 = verts(3,itrivert(2,itri)) *scale + dz
        v13 = verts(1,itrivert(3,itri)) *scale + dx
        v23 = verts(2,itrivert(3,itri)) *scale + dy
        v33 = verts(3,itrivert(3,itri)) *scale + dz
c        
        triangles(1,1,itri) = v11
        triangles(2,1,itri) = v21
        triangles(3,1,itri) = v31
        triangles(1,2,itri) = v12
        triangles(2,2,itri) = v22
        triangles(3,2,itri) = v32
        triangles(1,3,itri) = v13
        triangles(2,3,itri) = v23
        triangles(3,3,itri) = v33
        centroids(1,itri) = (v11+v12+v13)/3
        centroids(2,itri) = (v21+v22+v23)/3
        centroids(3,itri) = (v31+v32+v33)/3
        call triangle_norm(triangles(1,1,itri),trinorm(1,itri))
        call triangle_area(triangles(1,1,itri),triarea(itri))
c
        enddo
        return
        end
c
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
