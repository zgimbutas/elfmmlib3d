        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),du(3),rnorm(3),fvec(3),strain(3,3)
        dimension triangles(3,3,1000),trinorm(3,1000),sigma(3,1000)
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        rlam = 3.45d0
        rmu = 2.3d0
        rnu = rlam/2/(rlam+rmu)

        done=1 
        pi=4*atan(done)

        alpha = (rlam+rmu)/(rlam+2*rmu)

        call prin2('=== Point source (Okada) ===*',fvec,0)

        x=.2d0
        y=.3d0
        z=0

        depth = 1
        dip = 1
        phi = dip/180.0d0*pi

        pot1 = 1
        pot2 = 2
        pot3 = 3
        pot4 = 0

        ntimes=100000
        t1=second()
C$        t1=omp_get_wtime()
        do i=1,ntimes
        call DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,  
     $     UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET) 
        enddo
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('dc3d0: speed, particles/sec=*',ntimes/(t2-t1),1)

        fvec(1)=ux
        fvec(2)=uy
        fvec(3)=uz
        strain(1,1)=uxx
        strain(1,2)=(uxy+uyx)/2
        strain(1,3)=(uxz+uzx)/2
        strain(2,1)=(uyx+uxy)/2
        strain(2,2)=uyy
        strain(2,3)=(uyz+uzy)/2
        strain(3,1)=(uzx+uxz)/2
        strain(3,2)=(uzy+uyz)/2
        strain(3,3)=uzz
        call prin2('after dc3d0, fvec=*',fvec,3)
        call prin2('after dc3d0, strain=*',strain,3*3)

ccc        stop
        call prin2(' *',fvec,0)
        call prin2('=== Point source (Mindlin pieces) ===*',fvec,0)

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

        call prin2('du=*',du,3)
        call prin2('rnorm=*',rnorm,3)

        target(1)=x
        target(2)=y
        target(3)=z

        ifstrain = 1

        ntimes=100000
        t1=second()
C$        t1=omp_get_wtime()
        do i=1,ntimes
        call green3elth_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)
        enddo
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('green3elth_eval: speed, particles/sec=*',
     $     ntimes/(t2-t1),1)

        do i =1,3
        fvec(i)=fvec(i)/(4*pi)
        enddo
        do i =1,3
        do j =1,3
        strain(i,j)=strain(i,j)/(4*pi)
        enddo
        enddo

        call prin2('after green3elth_eval, fvec=*',fvec,3)
        call prin2('after green3elth_eval, strain=*',strain,3*3)
        
ccc     stop
        call prin2(' *',fvec,0)
        call prin2('=== Rectangle (Okada) ===*',fvec,0)

        al1 = 0
        al2 = 1
        aw1 = 0
        aw2 = 1
        disl1 = pot1
        disl2 = pot2
        disl3 = pot3

        ntimes=100000
        t1=second()
C$        t1=omp_get_wtime()
        do i=1,ntimes
        call  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,                        
     $     AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,           
     $     UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
        enddo
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('dc3d: speed, rectangles/sec=*',ntimes/(t2-t1),1)

        fvec(1)=ux
        fvec(2)=uy
        fvec(3)=uz
        strain(1,1)=uxx
        strain(1,2)=(uxy+uyx)/2
        strain(1,3)=(uxz+uzx)/2
        strain(2,1)=(uyx+uxy)/2
        strain(2,2)=uyy
        strain(2,3)=(uyz+uzy)/2
        strain(3,1)=(uzx+uxz)/2
        strain(3,2)=(uzy+uyz)/2
        strain(3,3)=uzz
        call prin2('after dc3d, fvec=*',fvec,3)
        call prin2('after dc3d, strain=*',strain,3*3)


ccc        stop
        call prin2(' *',fvec,0)
        call prin2('=== Triangles (adaptive integration) ===*',fvec,0)

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

        ntri = 2

        do i = 1,ntri
        sigma(1,i)=du(1)
        sigma(2,i)=du(2)
        sigma(3,i)=du(3)
        enddo

        call triangle_norm(triangles(1,1,1),trinorm(1,1))
        call triangle_norm(triangles(1,1,2),trinorm(1,2))
ccc        call prin2('triangles=*',triangles,3*3*ntri)
ccc        call prin2('trinorm=*',trinorm,3*ntri)

        ntimes=100
        t1=second()
C$        t1=omp_get_wtime()
        do i=1,ntimes
        call elthb3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,fvec,strain,numfunev)
        enddo
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('elthb3triaadap: speed, triangles/sec=*',
     $     ntimes*ntri/(t2-t1),1)
        call prinf('numfunev=*',numfunev,1)

        do i = 1,3
        fvec(i)=fvec(i)/(4*pi)
        enddo
        do i = 1,3
        do j = 1,3
        strain(i,j)=strain(i,j)/(4*pi)
        enddo
        enddo

        call prin2('after elthb3triaadap, fvec=*',fvec,3)
        call prin2('after elthb3triaadap, strain=*',strain,3*3)

        stop
        end
c
c
c
c
c
        subroutine elthb3triaadap
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
        external fun3elthb_eval
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
        
        nq = 6
        nfuns = 12
c       
        eps=1e-6
        call tria3adam(ier,vert1,vert2,vert3,fun3elthb_eval,nfuns,
     1      target,par,nq,eps,rints,maxrec,numfunev,w)
ccc        call prinf('ier=*',ier,1)
ccc        call prinf('numfunev=*',numfunev,1)

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
        subroutine fun3elthb_eval(x,y,z,target,par,f)
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
        call green3elth_eval
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

