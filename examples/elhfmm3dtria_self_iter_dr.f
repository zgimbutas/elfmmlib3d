c
c       This is a simple driver to test the layer potential FMM routines
c       in R^3 using the half-space Green's functions for general values 
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
        dimension ptfrc_direct(3,nmax),strain_direct(3,3,nmax)
        dimension ptfrc_directtarg(3,nmax),strain_directtarg(3,3,nmax)
c
        dimension ptfrc(3,nmax),strain(3,3,nmax)
        dimension ptfrctarg(3,nmax),straintarg(3,3,nmax)
c
        dimension ptfrc0(3),strain0(3,3),stress0(3,3),tract0(3)
c
        dimension ptfrc1(3,nmax),strain1(3,3,nmax)
        dimension ptfrc2(3,nmax),strain2(3,3,nmax)
c
        real *8, allocatable :: wsave_direct(:)
        real *8, allocatable :: wsave_himage(:)
c
c       SET ALL PARAMETERS
c
c       wsave_direct will be used to store local interactions
c       for the direct arrival from the source and 
c       wsave_himage will be used to store local interactions
c       for the image conrtibutions.
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
c       define piecewise constant densities
c
        do i=1,ntri
           sigma_sl(1,i)=1
           sigma_sl(2,i)=2
           sigma_sl(3,i)=3
           sigma_dl(1,i)=1
           sigma_dl(2,i)=2
           sigma_dl(3,i)=3
        enddo        
c
c       set targets on a regular grid near interface z=0
c
        ngrid=100
        kk=0
        do i=1,ngrid
        do j=1,ngrid
           kk=kk+1
           target(1,kk)=(i-1)
           target(2,kk)=(j-1)
           target(3,kk)=-1d-12
        enddo
        enddo
        ntargs=ngrid*ngrid
c
        ntargs=0
c
        call prinf('ntargs=*',ntargs,1)
ccc        call prin2('target=*',target,3*ntargs)
c
c
c     Check all sources and targets
c     For half space elastostatic Green's function, z MUST BE negative
c
        do i=1,ntri
           if( centroids(3,i) .gt. 0 ) then
              write(*,*) 'centroids(3,i) >=0 at source location',i
              call prin2('centroids=*',centroids(1,i),3)
              stop
           endif
        enddo
        do i=1,ntargs
           if( target(3,i) .gt. 0 ) then
              write(*,*) 'target(3,i) >= 0 at target location',i
              call prin2('target=*',target(1,i),3)
              stop
           endif
        enddo
c
        call prinf('==== targets ====*',i,0)
c       
c       since we are testing by direct calculation, only compute
c       direct calculation at modest number (m) of targets, set here.
c
ccc        m=min(ntargs,20)
ccc        call prinf('m=*',m,1)
c
c       ... m on sources only
c
        m=min(ntri,20)
        call prinf('m=*',m,1)
c
        rlam=3.45d0
        rmu=2.3d0
c
c       turn on single and/or double layer potential
c
        ifsingle=0
        ifdouble=1
c
c       set whether displacement and/or strain to be compute on surface
c       and at target locations.
c
        ifptfrc=0
        ifstrain=1
        ifptfrctarg=0
        ifstraintarg=0
c
c       ... evaluate via direct quadrature on triangles
c
        do i=1,ntri
        do j=1,3
           ptfrc2(j,i)=0
        enddo
        do j=1,3
        do k=1,3
           strain2(j,k,i)=0
        enddo
        enddo
        enddo
c
        do i=1,ntargs
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

        call test_elh3dtriadirecttarg(M,
     $     RLAM,RMU,TRIANGLES,TRINORM,NTRI,CENTROIDS,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc2,ifstrain,strain2,M,
     $     target,ifptfrctarg,PTFRC1,
     $     ifstraintarg,STRAIN1)

        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after elhfmm3triadirect, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (m)/(t2-t1),1)
        call prin2('after estimated time for direct=*',
     $     ntri/dble(m)*(t2-t1),1)
c
c-------------------------------------------------------------------
c       DIRECT ARRIVAL
c-------------------------------------------------------------------
c
c       ... evaluate direct arrival via FMM on triangles 
c
        iprec=2
c
c       One parameter to the FMM call is icomp_type.
c       Setting icomp_type = 4 estimates the amount of workspace
c       needed for the storage of local interactions.
c       It dones't actually do the FMM calculation (except for tree
c       generation.)
c
        icomp_type = 4
        t1=second()
C$        t1=omp_get_wtime()
        call elfmm3dtriatargiter
     $     (ier,iprec,RLAM,RMU,TRIANGLES,TRINORM,NTRI,centroids,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc_direct,ifstrain,strain_direct,
     $     ntargs,target,ifptfrctarg,PTFRC_directtarg,
     $     ifstraintarg,STRAIN_directtarg,
     $     icomp_type,wsave_direct,lwsave_direct,lused)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('icomp_type=4, lused(k)=*',lused/1000,1)
c       
c
c
        lwsave_direct = lused
c
c       allocated space for local interactions
c
        allocate( wsave_direct(lwsave_direct) )
c
c       When icomp_type = 1, the FMM is carried out without local 
c       storage.
c       When icomp_type = 2, the FMM is carried out AND it stores
c       local  interactions in wsave_direct.
c       storage.
c       When icomp_type = 3, the FMM is carried USING the local 
c       interactions stored in wsave_direct. It assumes a previous
c       call with icomp_type=2 was successfully completed.       
c       It requires the same FMM tree, so must be called with the same
c       parameters as the previous call to the FMM.
c
        do icomp_type = 2,3
c
        t1=second()
C$        t1=omp_get_wtime()
        call elfmm3dtriatargiter
     $     (ier,iprec,RLAM,RMU,TRIANGLES,TRINORM,NTRI,centroids,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc_direct,ifstrain,strain_direct,
     $     ntargs,target,ifptfrctarg,PTFRC_directtarg,
     $     ifstraintarg,STRAIN_directtarg,
     $     icomp_type,wsave_direct,lwsave_direct,lused)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after elhfmm3triatarg, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (ntargs)/(t2-t1),1)
        call prin2('speed, interactions(k)/sec=*',
     $     (ntargs*ntri)/1000/(t2-t1),1)
c
ccc        call prin2('wsave_direct=*',wsave_direct,2*120)
c
        enddo
c
ccc        call prin2('ptfrc_direct=*',ptfrc_direct,3*ntri)
ccc        call prin2('strain_direct=*',strain_direct,3*3*ntri)
c
c-------------------------------------------------------------------
c       IMAGE CONTRIBUTIONS
c-------------------------------------------------------------------
c
c       ... evaluate image contributions via FMM on triangles 
c
        icomp_type = 4
        t1=second()
C$        t1=omp_get_wtime()
        call elhfmm3dtriatargiter
     $     (ier,iprec,RLAM,RMU,TRIANGLES,TRINORM,NTRI,centroids,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ntri,centroids,ifptfrc,PTFRC,ifstrain,STRAIN,
     $     icomp_type,wsave_himage,lwsave_himage,lused)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('icomp_type=4, lused(k)=*',lused/1000,1)
c
        lwsave_himage = lused
        allocate( wsave_himage(lwsave_himage) )
c
        do icomp_type = 2,3
c
        t1=second()
C$        t1=omp_get_wtime()
        call elhfmm3dtriatargiter
     $     (ier,iprec,RLAM,RMU,TRIANGLES,TRINORM,NTRI,centroids,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ntri,centroids,ifptfrc,PTFRC,ifstrain,STRAIN,
     $     icomp_type,wsave_himage,lwsave_himage,lused)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after elhfmm3triatarg, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (ntargs)/(t2-t1),1)
        call prin2('speed, interactions(k)/sec=*',
     $     (ntargs*ntri)/1000/(t2-t1),1)
c
ccc        call prin2('wsave_himage=*',wsave_himage,2*120)
c
        enddo
c
ccc        call prin2('ptfrc_image=*',ptfrc,3*ntri)
ccc        call prin2('strain_image=*',strain,3*3*ntri)
c
c
c     ... add direct arrival to image contributions for
c         displacement/strain
c
        if( ifptfrc .eq. 1 ) then
        do i=1,ntri
         ptfrc(1,i)=ptfrc(1,i)+ptfrc_direct(1,i)
         ptfrc(2,i)=ptfrc(2,i)+ptfrc_direct(2,i)
         ptfrc(3,i)=ptfrc(3,i)+ptfrc_direct(3,i)
        enddo
        endif

        if( ifstrain .eq. 1 ) then
        do i=1,ntri
         strain(1,1,i)=strain(1,1,i)+strain_direct(1,1,i)
         strain(2,1,i)=strain(2,1,i)+strain_direct(2,1,i)
         strain(3,1,i)=strain(3,1,i)+strain_direct(3,1,i)
         strain(1,2,i)=strain(1,2,i)+strain_direct(1,2,i)
         strain(2,2,i)=strain(2,2,i)+strain_direct(2,2,i)
         strain(3,2,i)=strain(3,2,i)+strain_direct(3,2,i)
         strain(1,3,i)=strain(1,3,i)+strain_direct(1,3,i)
         strain(2,3,i)=strain(2,3,i)+strain_direct(2,3,i)
         strain(3,3,i)=strain(3,3,i)+strain_direct(3,3,i)
        enddo
        endif

        if( ifptfrctarg .eq. 1 ) then
        do i=1,ntargs
         ptfrctarg(1,i)=ptfrctarg(1,i)+ptfrc_directtarg(1,i)
         ptfrctarg(2,i)=ptfrctarg(2,i)+ptfrc_directtarg(2,i)
         ptfrctarg(3,i)=ptfrctarg(3,i)+ptfrc_directtarg(3,i)
        enddo
        endif

        if( ifstraintarg .eq. 1 ) then
        do i=1,ntargs
         straintarg(1,1,i)=straintarg(1,1,i)+strain_directtarg(1,1,i)
         straintarg(2,1,i)=straintarg(2,1,i)+strain_directtarg(2,1,i)
         straintarg(3,1,i)=straintarg(3,1,i)+strain_directtarg(3,1,i)
         straintarg(1,2,i)=straintarg(1,2,i)+strain_directtarg(1,2,i)
         straintarg(2,2,i)=straintarg(2,2,i)+strain_directtarg(2,2,i)
         straintarg(3,2,i)=straintarg(3,2,i)+strain_directtarg(3,2,i)
         straintarg(1,3,i)=straintarg(1,3,i)+strain_directtarg(1,3,i)
         straintarg(2,3,i)=straintarg(2,3,i)+strain_directtarg(2,3,i)
         straintarg(3,3,i)=straintarg(3,3,i)+strain_directtarg(3,3,i)
        enddo
        endif
c
c
c
        if( ifptfrc .eq. 1 ) 
     $     call prin2('after triadirect, ptfrc=*',ptfrc2,3*m)
        if( ifstrain .eq. 1 ) 
     $     call prin2('after triadirect, strain=*',strain2,3*3*m)
c
        if( ifptfrc .eq. 1 )
     $     call prin2('after elhfmm3tria, ptfrc=*',
     $     ptfrc,3*m)
        if( ifstrain .eq. 1 )
     $     call prin2('after elhfmm3tria, strain=*',
     $     strain,3*3*m)
c
c
        if( ifptfrctarg .eq. 1 ) 
     $   call prin2('after triadirecttarg, ptfrctarg=*',ptfrc1,3*m)
        if( ifstraintarg .eq. 1 ) 
     $   call prin2('after triadirecttarg, straintarg=*',strain1,3*3*m)
c
        if( ifptfrctarg .eq. 1 )
     $     call prin2('after elhfmm3triatarg, ptfrctarg=*',
     $     ptfrctarg,3*m)
        if( ifstraintarg .eq. 1 )
     $     call prin2('after elhfmm3triatarg, straintarg=*',
     $     straintarg,3*3*m)
c
c
        if( ifptfrc .eq. 1 ) then
        call d3error(ptfrc2,ptfrc,3*m,a,r)
        call prin2('absolute error in ptfrc=*',a,1)
        call prin2('relative error in ptfrc=*',r,1)
        endif
c
        if( ifstrain .eq. 1 ) then
        call d3error(strain2,strain,3*3*m,a,r)
        call prin2('absolute error in strain=*',a,1)
        call prin2('relative error in strain=*',r,1)
        endif
c
        if( ifptfrctarg .eq. 1 ) then
        call d3error(ptfrc1,ptfrctarg,3*m,a,r)
        call prin2('absolute error in target ptfrc=*',a,1)
        call prin2('relative error in target ptfrc=*',r,1)
        endif
c
        if( ifstraintarg .eq. 1 ) then
        call d3error(strain1,straintarg,3*3*m,a,r)
        call prin2('absolute error in target strain=*',a,1)
        call prin2('relative error in target strain=*',r,1)
        endif
c
c
        if (ifstraintarg .eq. 1) then
        do i=1,m
        call strain2stress(rlam,rmu,strain1(1,1,i),stress0)
        call prin2('directly, stresstarg=*',stress0,9)
        call strain2stress(rlam,rmu,straintarg(1,1,i),stress0)
        call prin2('after fmm, stresstarg=*',stress0,9)
        enddo
c
        endif
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
        open (unit = ir,file='./square.a.tri')

        call atrireadchk(ir,verts,lv,nverts,itrivert,li,ntri,iergeom)
c
c
ccc        ntri=179
ccc        ntri=12
ccc        ntri=340
c
        if (ntri.gt.nmax) then
        iergeom = 4
        return
        endif
c    
c       create triangles from tri format data
c
c       scale and offset if desired.
c
        scale=1
        dx=0
        dy=0
        dz=-1.5d0 
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
c
c
c
        subroutine strain2stress(rlam,rmu,strain,stress)
        implicit real *8 (a-h,o-z)
        dimension strain(3,3),stress(3,3)
c
c       ... convert strain tensor into the stress tensor
c
        do i=1,3
        do j=1,3
        stress(i,j)=0
        enddo
        enddo
c
        stress(1,1)=rlam*(strain(1,1)+strain(2,2)+strain(3,3))
        stress(2,2)=rlam*(strain(1,1)+strain(2,2)+strain(3,3))
        stress(3,3)=rlam*(strain(1,1)+strain(2,2)+strain(3,3))
c
        do i=1,3
        do j=1,3
        stress(i,j)=stress(i,j)+rmu*(strain(i,j)+strain(j,i))
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
        subroutine test_elh3dtriadirecttarg(M,
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     NTARGET,target,
     $     ifptfrctarg,PTFRCtarg,ifstraintarg,STRAINtarg)
        implicit real *8 (a-h,o-z)
c
c
c       Elastostatic interactions in R^3: evaluate all pairwise triangle
c       interactions and interactions with targets.  In order to
c       accelerate tests, self strains and displacements are evaluated
c       for first m source locations only.
c
c       INPUT:
c
c       M - the number of source locations for testing
c       rlam,rmu - Lame parameters
c       triangle(3,3,nsource) - array of triangles in standard format
c       trinorm(3,nsource) - array of triangle normals
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifptfrc - displacement computation flag
c       ifstrain - strain computation flag
c       ifptfrctarg - target displacement computation flag
c       ifstraintarg - target strain computation flag
c
c       OUTPUT:
c
c       ptfrc(3,m) - displacement at source locations
c       strain(3,3,m) - strain at source locations
c       ptfrctarg(3,ntarget) - displacement at target locations
c       straintarg(3,3,ntarget) - strain at target locations
c
c
        dimension triangle(3,3,1),trinorm(3,1),source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension ptfrc(3,1),strain(3,3,1)
        dimension ptfrctarg(3,1),straintarg(3,3,1)
c
        dimension ptfrc0(3),strain0(3,3)
c
c
c     NOTE: In the present version, dipole vectors SIGMA_DV must be SET EQUAL
c     to the triangle normal for elt3d direct routines
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
        ione=1
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0,numfunev)
        do j=1,m
        do i=1,nsource
c
        if (ifsingle .eq. 1 ) then
        if( i .eq. j ) then
        call elust3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        else
        call elust3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
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
        if( i .eq. j ) then 
        call eltst3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        else
        call eltst3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
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
        enddo
c
c       ... image contribution
c
        if (ifsingle .eq. 1 ) then
        call eluh3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,trinorm,source(1,j),
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if( ifstrain .eq. 1 ) then
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

        if (ifdouble .eq. 1 ) then
        call elth3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,source(1,j),
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if( ifstrain .eq. 1 ) then
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

        enddo
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptfrc0,strain0,numfunev)
        do j=1,ntarget
        if (ifsingle .eq. 1 ) then
        call elust3triadirecttarg
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
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
        if (ifdouble .eq. 1) then
        call eltst3triadirecttarg
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
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
c       ... image contribution
c
        if (ifsingle .eq. 1 ) then
        call eluh3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,trinorm,target(1,j),
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)
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

        if (ifdouble .eq. 1 ) then
        call elth3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,target(1,j),
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)
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
