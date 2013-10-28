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
c    $Date: 2011-07-13 21:58:51 -0400 (Wed, 13 Jul 2011) $
c    $Revision: 2226 $
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of various (free space) elastostatic 
c       Green's functions in R^3.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine green3elt(rlam,rmu,xyz,du,rnorm,fvec)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec at the
c       location xyz due to the static double force (aka slip,
c       displacement jump) du at the origin with the orientation vector
c       in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       displacement can be in any direction.
c       
c       Formula (38) in T. Maruyama's paper "On the Force Equivalents 
c       of Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, Vol. 41 
c       (1963), pp. 467-486.
c
c       Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       xyz (real *8) - the target point in R^3
c       du (real *8) - the strength of the double force source
c       rnorm (real *8) - the orientation vector of the double force source
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c
c
        dimension xyz(3),du(3),rnorm(3),fvec(3),ftmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        ftmp(1)=-dd*dx
        ftmp(2)=-dd*dy
        ftmp(3)=-dd*dz
c       
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        ftmp(1)=ftmp(1)+dd*du(1)
        ftmp(2)=ftmp(2)+dd*du(2)
        ftmp(3)=ftmp(3)+dd*du(3)
c       
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        ftmp(1)=ftmp(1)+dd*rnorm(1)
        ftmp(2)=ftmp(2)+dd*rnorm(2)
        ftmp(3)=ftmp(3)+dd*rnorm(3)
c       
        c1=rmu/(rlam+2*rmu)
c
        fvec(1)=ftmp(1)*c1/cd**3
        fvec(2)=ftmp(2)*c1/cd**3
        fvec(3)=ftmp(3)*c1/cd**3
c
c       ... step 2
c
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
c
        c2=(rlam+rmu)/(rlam+2*rmu)
c
        fvec(1)=fvec(1)+3*c2*dd*dx/cd**5
        fvec(2)=fvec(2)+3*c2*dd*dy/cd**5
        fvec(3)=fvec(3)+3*c2*dd*dz/cd**5
c
        return
        end
c
c
c
c
c
        subroutine green3elt_strain(rlam,rmu,xyz,du,rnorm,fvec,strain)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the double force
c       (aka slip, displacement jump) du at the origin with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Formula (38) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       xyz (real *8) - the target point in R^3
c       du (real *8) - the strength of the double force source
c       rnorm (real *8) - the orientation vector of the double force source
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c       strain (real *8) - the strain at the target
c
c
        dimension xyz(3),du(3),rnorm(3),fvec(3),ftmp(3)
        dimension strain(3,3),dudx(3,3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        ftmp(1)=-dd*dx
        ftmp(2)=-dd*dy
        ftmp(3)=-dd*dz
c       
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        ftmp(1)=ftmp(1)+dd*du(1)
        ftmp(2)=ftmp(2)+dd*du(2)
        ftmp(3)=ftmp(3)+dd*du(3)
c       
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        ftmp(1)=ftmp(1)+dd*rnorm(1)
        ftmp(2)=ftmp(2)+dd*rnorm(2)
        ftmp(3)=ftmp(3)+dd*rnorm(3)
c       
        c1=rmu/(rlam+2*rmu)
c
        fvec(1)=ftmp(1)*c1/cd**3
        fvec(2)=ftmp(2)*c1/cd**3
        fvec(3)=ftmp(3)*c1/cd**3
c
        dudx(1,1)=ftmp(1)*c1*(-3*dx/cd**5)
        dudx(1,2)=ftmp(1)*c1*(-3*dy/cd**5)
        dudx(1,3)=ftmp(1)*c1*(-3*dz/cd**5)

        dudx(2,1)=ftmp(2)*c1*(-3*dx/cd**5)
        dudx(2,2)=ftmp(2)*c1*(-3*dy/cd**5)
        dudx(2,3)=ftmp(2)*c1*(-3*dz/cd**5)

        dudx(3,1)=ftmp(3)*c1*(-3*dx/cd**5)
        dudx(3,2)=ftmp(3)*c1*(-3*dy/cd**5)
        dudx(3,3)=ftmp(3)*c1*(-3*dz/cd**5)

        dudx(1,1)=dudx(1,1)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)) * c1/cd**3
        dudx(2,2)=dudx(2,2)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)) * c1/cd**3
        dudx(3,3)=dudx(3,3)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)) * c1/cd**3

        dudx(1,1)=dudx(1,1)+(du(1)*rnorm(1)) * c1/cd**3
        dudx(1,2)=dudx(1,2)+(du(1)*rnorm(2)) * c1/cd**3
        dudx(1,3)=dudx(1,3)+(du(1)*rnorm(3)) * c1/cd**3
        dudx(2,1)=dudx(2,1)+(du(2)*rnorm(1)) * c1/cd**3
        dudx(2,2)=dudx(2,2)+(du(2)*rnorm(2)) * c1/cd**3
        dudx(2,3)=dudx(2,3)+(du(2)*rnorm(3)) * c1/cd**3
        dudx(3,1)=dudx(3,1)+(du(3)*rnorm(1)) * c1/cd**3
        dudx(3,2)=dudx(3,2)+(du(3)*rnorm(2)) * c1/cd**3
        dudx(3,3)=dudx(3,3)+(du(3)*rnorm(3)) * c1/cd**3
c
        dudx(1,1)=dudx(1,1)+(du(1)*rnorm(1)) * c1/cd**3
        dudx(1,2)=dudx(1,2)+(du(2)*rnorm(1)) * c1/cd**3
        dudx(1,3)=dudx(1,3)+(du(3)*rnorm(1)) * c1/cd**3
        dudx(2,1)=dudx(2,1)+(du(1)*rnorm(2)) * c1/cd**3
        dudx(2,2)=dudx(2,2)+(du(2)*rnorm(2)) * c1/cd**3
        dudx(2,3)=dudx(2,3)+(du(3)*rnorm(2)) * c1/cd**3
        dudx(3,1)=dudx(3,1)+(du(1)*rnorm(3)) * c1/cd**3
        dudx(3,2)=dudx(3,2)+(du(2)*rnorm(3)) * c1/cd**3
        dudx(3,3)=dudx(3,3)+(du(3)*rnorm(3)) * c1/cd**3
c
c       ... step 2
c
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
c
        c2=(rlam+rmu)/(rlam+2*rmu)
c
        fvec(1)=fvec(1)+3*c2*dd*dx/cd**5
        fvec(2)=fvec(2)+3*c2*dd*dy/cd**5
        fvec(3)=fvec(3)+3*c2*dd*dz/cd**5
c
        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
c
        dudx(1,1)=dudx(1,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dx/cd**5
        dudx(1,2)=dudx(1,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dx/cd**5
        dudx(1,3)=dudx(1,3)+3*c2*(du(3)*d2+rnorm(3)*d1)*dx/cd**5
        dudx(2,1)=dudx(2,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dy/cd**5
        dudx(2,2)=dudx(2,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dy/cd**5
        dudx(2,3)=dudx(2,3)+3*c2*(du(3)*d2+rnorm(3)*d1)*dy/cd**5
        dudx(3,1)=dudx(3,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dz/cd**5
        dudx(3,2)=dudx(3,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dz/cd**5
        dudx(3,3)=dudx(3,3)+3*c2*(du(3)*d2+rnorm(3)*d1)*dz/cd**5
c
        dudx(1,1)=dudx(1,1)+3*c2*dd*(1/cd**5-5*dx*dx/cd**7)
        dudx(1,2)=dudx(1,2)+3*c2*dd*(       -5*dx*dy/cd**7)
        dudx(1,3)=dudx(1,3)+3*c2*dd*(       -5*dx*dz/cd**7)
        dudx(2,1)=dudx(2,1)+3*c2*dd*(       -5*dy*dx/cd**7)
        dudx(2,2)=dudx(2,2)+3*c2*dd*(1/cd**5-5*dy*dy/cd**7)
        dudx(2,3)=dudx(2,3)+3*c2*dd*(       -5*dy*dz/cd**7)
        dudx(3,1)=dudx(3,1)+3*c2*dd*(       -5*dz*dx/cd**7)
        dudx(3,2)=dudx(3,2)+3*c2*dd*(       -5*dz*dy/cd**7)
        dudx(3,3)=dudx(3,3)+3*c2*dd*(1/cd**5-5*dz*dz/cd**7)
c
c       ... evaluate strain
c       
        do i=1,3
        do j=1,3
        strain(i,j)=(dudx(i,j)+dudx(j,i))/2
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
        subroutine green3elt_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location target due to the double force
c       (aka slip, displacement jump) du at source with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Formula (38) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       source (real *8) - the source point in R^3
c       du (real *8) - the strength of the double force source
c       rnorm (real *8) - the orientation vector of the double force source
c       target (real *8) - the target point in R^3
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c       strain (real *8) - the strain at the target
c
c
        dimension source(3),target(3),du(3),rnorm(3)
        dimension fvec(3),ftmp(3)
        dimension strain(3,3),dudx(3,3)
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
c
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
        c1=rmu/(rlam+2*rmu)
c
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        ftmp(1)=-dd*dx
        ftmp(2)=-dd*dy
        ftmp(3)=-dd*dz
c       
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        ftmp(1)=ftmp(1)+dd*du(1)
        ftmp(2)=ftmp(2)+dd*du(2)
        ftmp(3)=ftmp(3)+dd*du(3)
c       
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        ftmp(1)=ftmp(1)+dd*rnorm(1)
        ftmp(2)=ftmp(2)+dd*rnorm(2)
        ftmp(3)=ftmp(3)+dd*rnorm(3)
c       
        fvec(1)=ftmp(1)*c1/cd**3
        fvec(2)=ftmp(2)*c1/cd**3
        fvec(3)=ftmp(3)*c1/cd**3
c
        if( ifstrain .eq. 1 ) then
c
        dudx(1,1)=ftmp(1)*c1*(-3*dx/cd**5)
        dudx(1,2)=ftmp(1)*c1*(-3*dy/cd**5)
        dudx(1,3)=ftmp(1)*c1*(-3*dz/cd**5)

        dudx(2,1)=ftmp(2)*c1*(-3*dx/cd**5)
        dudx(2,2)=ftmp(2)*c1*(-3*dy/cd**5)
        dudx(2,3)=ftmp(2)*c1*(-3*dz/cd**5)

        dudx(3,1)=ftmp(3)*c1*(-3*dx/cd**5)
        dudx(3,2)=ftmp(3)*c1*(-3*dy/cd**5)
        dudx(3,3)=ftmp(3)*c1*(-3*dz/cd**5)

        dudx(1,1)=dudx(1,1)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3))*c1/cd**3
        dudx(2,2)=dudx(2,2)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3))*c1/cd**3
        dudx(3,3)=dudx(3,3)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3))*c1/cd**3

        dudx(1,1)=dudx(1,1)+(du(1)*rnorm(1)) * c1/cd**3
        dudx(1,2)=dudx(1,2)+(du(1)*rnorm(2)) * c1/cd**3
        dudx(1,3)=dudx(1,3)+(du(1)*rnorm(3)) * c1/cd**3
        dudx(2,1)=dudx(2,1)+(du(2)*rnorm(1)) * c1/cd**3
        dudx(2,2)=dudx(2,2)+(du(2)*rnorm(2)) * c1/cd**3
        dudx(2,3)=dudx(2,3)+(du(2)*rnorm(3)) * c1/cd**3
        dudx(3,1)=dudx(3,1)+(du(3)*rnorm(1)) * c1/cd**3
        dudx(3,2)=dudx(3,2)+(du(3)*rnorm(2)) * c1/cd**3
        dudx(3,3)=dudx(3,3)+(du(3)*rnorm(3)) * c1/cd**3
c
        dudx(1,1)=dudx(1,1)+(du(1)*rnorm(1)) * c1/cd**3
        dudx(1,2)=dudx(1,2)+(du(2)*rnorm(1)) * c1/cd**3
        dudx(1,3)=dudx(1,3)+(du(3)*rnorm(1)) * c1/cd**3
        dudx(2,1)=dudx(2,1)+(du(1)*rnorm(2)) * c1/cd**3
        dudx(2,2)=dudx(2,2)+(du(2)*rnorm(2)) * c1/cd**3
        dudx(2,3)=dudx(2,3)+(du(3)*rnorm(2)) * c1/cd**3
        dudx(3,1)=dudx(3,1)+(du(1)*rnorm(3)) * c1/cd**3
        dudx(3,2)=dudx(3,2)+(du(2)*rnorm(3)) * c1/cd**3
        dudx(3,3)=dudx(3,3)+(du(3)*rnorm(3)) * c1/cd**3
c
        endif
c
c       ... step 2
c
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
c
        c2=(rlam+rmu)/(rlam+2*rmu)
c
        fvec(1)=fvec(1)+3*c2*dd*dx/cd**5
        fvec(2)=fvec(2)+3*c2*dd*dy/cd**5
        fvec(3)=fvec(3)+3*c2*dd*dz/cd**5
c
        if( ifstrain .eq. 1 ) then
c
        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
c
        dudx(1,1)=dudx(1,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dx/cd**5
        dudx(1,2)=dudx(1,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dx/cd**5
        dudx(1,3)=dudx(1,3)+3*c2*(du(3)*d2+rnorm(3)*d1)*dx/cd**5
        dudx(2,1)=dudx(2,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dy/cd**5
        dudx(2,2)=dudx(2,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dy/cd**5
        dudx(2,3)=dudx(2,3)+3*c2*(du(3)*d2+rnorm(3)*d1)*dy/cd**5
        dudx(3,1)=dudx(3,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dz/cd**5
        dudx(3,2)=dudx(3,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dz/cd**5
        dudx(3,3)=dudx(3,3)+3*c2*(du(3)*d2+rnorm(3)*d1)*dz/cd**5
c
        dudx(1,1)=dudx(1,1)+3*c2*dd*(1/cd**5-5*dx*dx/cd**7)
        dudx(1,2)=dudx(1,2)+3*c2*dd*(       -5*dx*dy/cd**7)
        dudx(1,3)=dudx(1,3)+3*c2*dd*(       -5*dx*dz/cd**7)
        dudx(2,1)=dudx(2,1)+3*c2*dd*(       -5*dy*dx/cd**7)
        dudx(2,2)=dudx(2,2)+3*c2*dd*(1/cd**5-5*dy*dy/cd**7)
        dudx(2,3)=dudx(2,3)+3*c2*dd*(       -5*dy*dz/cd**7)
        dudx(3,1)=dudx(3,1)+3*c2*dd*(       -5*dz*dx/cd**7)
        dudx(3,2)=dudx(3,2)+3*c2*dd*(       -5*dz*dy/cd**7)
        dudx(3,3)=dudx(3,3)+3*c2*dd*(1/cd**5-5*dz*dz/cd**7)
c
        endif
c
c
        if( ifstrain .eq. 1 ) then
c
c       ... evaluate strain
c       
        do i=1,3
        do j=1,3
        strain(i,j)=(dudx(i,j)+dudx(j,i))/2
        enddo
        enddo
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
        subroutine green3elu(rlam,rmu,xyz,df,fvec)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic single layer Green's function: Kelvin solution
c
c       This subroutine evaluates the displacement vector fvec at the
c       location xyz due to the static single force (aka traction) df
c       at the origin.
c
c       Formula (39) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       source (real *8 ) - the source point in R^3
c       df (real *8) - the strength of the single force source 
c       target (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c
c
        dimension xyz(3),dtau(3,3),rnorm(3),fvec(3),df(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
        c2=(rlam+rmu)/(rlam+2*rmu)
ccc        c1=(rlam+3*rmu)/(rlam+2*rmu)
        c1 = 2.0d0-c2
c
        fvec(1)=df(1)*c1/cd
        fvec(2)=df(2)*c1/cd
        fvec(3)=df(3)*c1/cd
c
c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
c
        fvec(1)=fvec(1)+c2*dd*dx/cd**3
        fvec(2)=fvec(2)+c2*dd*dy/cd**3
        fvec(3)=fvec(3)+c2*dd*dz/cd**3
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(2*rmu)
        fvec(2)=fvec(2)/(2*rmu)
        fvec(3)=fvec(3)/(2*rmu)
c
        return
        end
c
c
c
c
c 
        subroutine green3elu_strain(rlam,rmu,xyz,df,fvec,strain)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic single layer Green's function: Kelvin solution
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the static
c       single force (aka traction) df at the origin.
c
c       Formula (39) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       xyz (real *8 ) - the target point in R^3
c       df (real *8) - the strength of the single force source 
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c       strain (real *8) - the strain at the target
c
c
        dimension xyz(3),fvec(3),df(3)
        dimension strain(3,3),dudx(3,3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
ccc        c1=(rlam+3*rmu)/(rlam+2*rmu)
        c2=(rlam+rmu)/(rlam+2*rmu)
        c1=2.0d0-c2
ccc        call prin2(' c1 is *',c1,1)
ccc        call prin2(' c2 is *',c2,1)
c
        fvec(1)=df(1)*c1/cd
        fvec(2)=df(2)*c1/cd
        fvec(3)=df(3)*c1/cd
c
        dudx(1,1)=-df(1)*c1*dx/cd**3
        dudx(1,2)=-df(1)*c1*dy/cd**3
        dudx(1,3)=-df(1)*c1*dz/cd**3
        dudx(2,1)=-df(2)*c1*dx/cd**3
        dudx(2,2)=-df(2)*c1*dy/cd**3
        dudx(2,3)=-df(2)*c1*dz/cd**3
        dudx(3,1)=-df(3)*c1*dx/cd**3
        dudx(3,2)=-df(3)*c1*dy/cd**3
        dudx(3,3)=-df(3)*c1*dz/cd**3
c
c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
c
        fvec(1)=fvec(1)+c2*dd*dx/cd**3
        fvec(2)=fvec(2)+c2*dd*dy/cd**3
        fvec(3)=fvec(3)+c2*dd*dz/cd**3
c
        dudx(1,1)=dudx(1,1)+c2*df(1)*(dx/cd**3)
        dudx(1,2)=dudx(1,2)+c2*df(2)*(dx/cd**3)
        dudx(1,3)=dudx(1,3)+c2*df(3)*(dx/cd**3)
        dudx(2,1)=dudx(2,1)+c2*df(1)*(dy/cd**3)
        dudx(2,2)=dudx(2,2)+c2*df(2)*(dy/cd**3)
        dudx(2,3)=dudx(2,3)+c2*df(3)*(dy/cd**3)
        dudx(3,1)=dudx(3,1)+c2*df(1)*(dz/cd**3)
        dudx(3,2)=dudx(3,2)+c2*df(2)*(dz/cd**3)
        dudx(3,3)=dudx(3,3)+c2*df(3)*(dz/cd**3)
c
        dudx(1,1)=dudx(1,1)+c2*dd*(1/cd**3-3*dx*dx/cd**5)
        dudx(1,2)=dudx(1,2)+c2*dd*(       -3*dx*dy/cd**5)
        dudx(1,3)=dudx(1,3)+c2*dd*(       -3*dx*dz/cd**5)
        dudx(2,1)=dudx(2,1)+c2*dd*(       -3*dy*dx/cd**5)
        dudx(2,2)=dudx(2,2)+c2*dd*(1/cd**3-3*dy*dy/cd**5)
        dudx(2,3)=dudx(2,3)+c2*dd*(       -3*dy*dz/cd**5)
        dudx(3,1)=dudx(3,1)+c2*dd*(       -3*dz*dx/cd**5)
        dudx(3,2)=dudx(3,2)+c2*dd*(       -3*dz*dy/cd**5)
        dudx(3,3)=dudx(3,3)+c2*dd*(1/cd**3-3*dz*dz/cd**5)
c
c       ... evaluate strain
c       
        do i=1,3
        do j=1,3
        strain(i,j)=(dudx(i,j)+dudx(j,i))/2
        enddo
        enddo
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(2*rmu)
        fvec(2)=fvec(2)/(2*rmu)
        fvec(3)=fvec(3)/(2*rmu)
c
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
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
c
        subroutine green3elu_eval(rlam,rmu,source,df,target,
     $     fvec,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic single layer Green's function: Kelvin solution
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location target due to the static
c       single force (aka traction) df at source.
c
c       Formula (39) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       xyz (real *8 ) - the target point in R^3
c       df (real *8) - the strength of the single force source 
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c       strain (real *8) - the strain at the target
c
c
        dimension source(3),target(3),df(3),fvec(3)
        dimension strain(3,3),dudx(3,3)
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
c
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
ccc        c1=(rlam+3*rmu)/(rlam+2*rmu)
        c2=(rlam+rmu)/(rlam+2*rmu)
        c1=2.0d0-c2
ccc        call prin2(' c1 is *',c1,1)
ccc        call prin2(' c2 is *',c2,1)
c
        fvec(1)=df(1)*c1/cd
        fvec(2)=df(2)*c1/cd
        fvec(3)=df(3)*c1/cd
c
        if( ifstrain .eq. 1 ) then
c
        dudx(1,1)=-df(1)*c1*dx/cd**3
        dudx(1,2)=-df(1)*c1*dy/cd**3
        dudx(1,3)=-df(1)*c1*dz/cd**3
        dudx(2,1)=-df(2)*c1*dx/cd**3
        dudx(2,2)=-df(2)*c1*dy/cd**3
        dudx(2,3)=-df(2)*c1*dz/cd**3
        dudx(3,1)=-df(3)*c1*dx/cd**3
        dudx(3,2)=-df(3)*c1*dy/cd**3
        dudx(3,3)=-df(3)*c1*dz/cd**3
c
        endif
c
c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
c
        fvec(1)=fvec(1)+c2*dd*dx/cd**3
        fvec(2)=fvec(2)+c2*dd*dy/cd**3
        fvec(3)=fvec(3)+c2*dd*dz/cd**3
c
        if( ifstrain .eq. 1 ) then
c
        dudx(1,1)=dudx(1,1)+c2*df(1)*(dx/cd**3)
        dudx(1,2)=dudx(1,2)+c2*df(2)*(dx/cd**3)
        dudx(1,3)=dudx(1,3)+c2*df(3)*(dx/cd**3)
        dudx(2,1)=dudx(2,1)+c2*df(1)*(dy/cd**3)
        dudx(2,2)=dudx(2,2)+c2*df(2)*(dy/cd**3)
        dudx(2,3)=dudx(2,3)+c2*df(3)*(dy/cd**3)
        dudx(3,1)=dudx(3,1)+c2*df(1)*(dz/cd**3)
        dudx(3,2)=dudx(3,2)+c2*df(2)*(dz/cd**3)
        dudx(3,3)=dudx(3,3)+c2*df(3)*(dz/cd**3)
c
        dudx(1,1)=dudx(1,1)+c2*dd*(1/cd**3-3*dx*dx/cd**5)
        dudx(1,2)=dudx(1,2)+c2*dd*(       -3*dx*dy/cd**5)
        dudx(1,3)=dudx(1,3)+c2*dd*(       -3*dx*dz/cd**5)
        dudx(2,1)=dudx(2,1)+c2*dd*(       -3*dy*dx/cd**5)
        dudx(2,2)=dudx(2,2)+c2*dd*(1/cd**3-3*dy*dy/cd**5)
        dudx(2,3)=dudx(2,3)+c2*dd*(       -3*dy*dz/cd**5)
        dudx(3,1)=dudx(3,1)+c2*dd*(       -3*dz*dx/cd**5)
        dudx(3,2)=dudx(3,2)+c2*dd*(       -3*dz*dy/cd**5)
        dudx(3,3)=dudx(3,3)+c2*dd*(1/cd**3-3*dz*dz/cd**5)
c
        endif
c
        if( ifstrain .eq. 1 ) then
c
c       ... evaluate strain
c       
        do i=1,3
        do j=1,3
        strain(i,j)=(dudx(i,j)+dudx(j,i))/2
        enddo
        enddo
c
        endif
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(2*rmu)
        fvec(2)=fvec(2)/(2*rmu)
        fvec(3)=fvec(3)/(2*rmu)
c
        if( ifstrain .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
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
        subroutine green3elu_strain_image(rlam,rmu,source,df,target,
     1             fvec,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic single layer Green's function: image A
c
c       ------------------------
c       In the Mindlin/Okada solution, there is an image
c       calculation that looks like the field due to the original
c       source at the image of the target. fvec is computed as usual
c       and derivatives w.r.t x,y are the same, but every z derivative
c       needs a sign flip before assembling the strain matrix.
c       ------------------------
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the static
c       single force (aka traction) df at the origin.
c
c       Formula (39) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       xyz (real *8 ) - the target point in R^3
c       df (real *8) - the strength of the single force source 
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c       strain (real *8) - the strain at the target
c
c
        dimension source(3),target(3)
        dimension xyz(3),dtau(3,3),rnorm(3),fvec(3),df(3)
        dimension strain(3,3),dudx(3,3)
c
        xyz(1) = target(1)-source(1)
        xyz(2) = target(2)-source(2)
        xyz(3) = target(3)-source(3)
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
ccc        c1=(rlam+3*rmu)/(rlam+2*rmu)
        c2=(rlam+rmu)/(rlam+2*rmu)
        c1=2.0d0-c2
ccc        call prin2(' c1 is *',c1,1)
ccc        call prin2(' c2 is *',c2,1)
c
        fvec(1)=df(1)*c1/cd
        fvec(2)=df(2)*c1/cd
        fvec(3)=df(3)*c1/cd
c
      if (ifstrain.eq.1) then
         dudx(1,1)=-df(1)*c1*dx/cd**3
         dudx(1,2)=-df(1)*c1*dy/cd**3
         dudx(1,3)=+df(1)*c1*dz/cd**3
         dudx(2,1)=-df(2)*c1*dx/cd**3
         dudx(2,2)=-df(2)*c1*dy/cd**3
         dudx(2,3)=+df(2)*c1*dz/cd**3
         dudx(3,1)=-df(3)*c1*dx/cd**3
         dudx(3,2)=-df(3)*c1*dy/cd**3
         dudx(3,3)=+df(3)*c1*dz/cd**3
      endif
c
c       ... step 2
c
      dd=df(1)*dx+df(2)*dy+df(3)*dz
c
      fvec(1)=fvec(1)+c2*dd*dx/cd**3
      fvec(2)=fvec(2)+c2*dd*dy/cd**3
      fvec(3)=fvec(3)+c2*dd*dz/cd**3
c
      if (ifstrain.eq.1) then
         dudx(1,1)=dudx(1,1)+c2*df(1)*(dx/cd**3)
         dudx(1,2)=dudx(1,2)+c2*df(2)*(dx/cd**3)
         dudx(1,3)=dudx(1,3)-c2*df(3)*(dx/cd**3)
         dudx(2,1)=dudx(2,1)+c2*df(1)*(dy/cd**3)
         dudx(2,2)=dudx(2,2)+c2*df(2)*(dy/cd**3)
         dudx(2,3)=dudx(2,3)-c2*df(3)*(dy/cd**3)
         dudx(3,1)=dudx(3,1)+c2*df(1)*(dz/cd**3)
         dudx(3,2)=dudx(3,2)+c2*df(2)*(dz/cd**3)
         dudx(3,3)=dudx(3,3)-c2*df(3)*(dz/cd**3)
c
         dudx(1,1)=dudx(1,1)+c2*dd*(1/cd**3-3*dx*dx/cd**5)
         dudx(1,2)=dudx(1,2)+c2*dd*(       -3*dx*dy/cd**5)
         dudx(1,3)=dudx(1,3)-c2*dd*(       -3*dx*dz/cd**5)
         dudx(2,1)=dudx(2,1)+c2*dd*(       -3*dy*dx/cd**5)
         dudx(2,2)=dudx(2,2)+c2*dd*(1/cd**3-3*dy*dy/cd**5)
         dudx(2,3)=dudx(2,3)-c2*dd*(       -3*dy*dz/cd**5)
         dudx(3,1)=dudx(3,1)+c2*dd*(       -3*dz*dx/cd**5)
         dudx(3,2)=dudx(3,2)+c2*dd*(       -3*dz*dy/cd**5)
         dudx(3,3)=dudx(3,3)-c2*dd*(1/cd**3-3*dz*dz/cd**5)
c
c       ... evaluate strain
c       
         do i=1,3
         do j=1,3
         strain(i,j)=(dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
      endif
c
c       ... scale by 1/(2*rmu)
c
      fvec(1)=fvec(1)/(2*rmu)
      fvec(2)=fvec(2)/(2*rmu)
      fvec(3)=fvec(3)/(2*rmu)
c
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
         strain(i,j)=strain(i,j)/(2*rmu)
         enddo
         enddo
      endif
c
      return
      end
c
c
c
c
        subroutine green3elt_strain_image(rlam,rmu,alpha,source,du,
     1             rnorm,target,fvec,ifstrain,strain)
        implicit real *8 (a-h,o-z)
c
c       Elastostatic double layer Green's function: image A
c
c       ------------------------
c       In the Mindlin/Okada solution, there is an image
c       calculation that looks like the field due to the original
c       source at the image of the target. fvec is computed as usual
c       and derivatives w.r.t x,y are the same, but every z derivative
c       needs a sign flip before assembling the strain matrix.
c       ------------------------
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the double force
c       (aka slip, displacement jump) du at the origin with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Formula (38) in T. Maruyama's paper "On the Force Equivalents of
c       Dynamic Elastic Dislocations with Reference to Earthquake 
c       Mechanism," Bull. of the Earthquake Research Institute, 
c       Vol. 41 (1963), pp. 467-486.
c
c
c          Input parameters:
c
c       rlam, rmu (real *8) - Lame parameters
c       xyz (real *8) - the target point in R^3
c       du (real *8) - the strength of the double force source
c       rnorm (real *8) - the orientation vector of the double force source
c
c          Output parameters:
c
c       fvec (real *8) - the displacement at the target
c       strain (real *8) - the strain at the target
c
c
        dimension source(3),target(3)
        dimension xyz(3),du(3),rnorm(3),fvec(3),ftmp(3)
        dimension strain(3,3),dudx(3,3)
c
        xyz(1) = target(1)-source(1)
        xyz(2) = target(2)-source(2)
        xyz(3) = target(3)-source(3)
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
	alphaim = (2-alpha)
        c1=1-alphaim
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        ftmp(1)=-dd*dx
        ftmp(2)=-dd*dy
        ftmp(3)=-dd*dz
c       
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        ftmp(1)=ftmp(1)+dd*du(1)
        ftmp(2)=ftmp(2)+dd*du(2)
        ftmp(3)=ftmp(3)+dd*du(3)
c       
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        ftmp(1)=ftmp(1)+dd*rnorm(1)
        ftmp(2)=ftmp(2)+dd*rnorm(2)
        ftmp(3)=ftmp(3)+dd*rnorm(3)
c       
        fvec(1)=ftmp(1)*c1/cd**3
        fvec(2)=ftmp(2)*c1/cd**3
        fvec(3)=ftmp(3)*c1/cd**3
c
ccc        call prin2(' inside image routine fvec is *',fvec,3)
        do i = 1,3
        do j = 1,3
           dudx(i,j) = 0.0d0
	enddo
	enddo
        if( ifstrain .eq. 1 ) then
c
        dudx(1,1)=ftmp(1)*c1*(-3*dx/cd**5)
        dudx(1,2)=ftmp(1)*c1*(-3*dy/cd**5)
        dudx(1,3)=-ftmp(1)*c1*(-3*dz/cd**5)

        dudx(2,1)=ftmp(2)*c1*(-3*dx/cd**5)
        dudx(2,2)=ftmp(2)*c1*(-3*dy/cd**5)
        dudx(2,3)=-ftmp(2)*c1*(-3*dz/cd**5)

        dudx(3,1)=ftmp(3)*c1*(-3*dx/cd**5)
        dudx(3,2)=ftmp(3)*c1*(-3*dy/cd**5)
        dudx(3,3)=-ftmp(3)*c1*(-3*dz/cd**5)

        dudx(1,1)=dudx(1,1)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3))*c1/cd**3
        dudx(2,2)=dudx(2,2)-
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3))*c1/cd**3
        dudx(3,3)=dudx(3,3)+
     $     (du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3))*c1/cd**3

        dudx(1,1)=dudx(1,1)+(du(1)*rnorm(1)) * c1/cd**3
        dudx(1,2)=dudx(1,2)+(du(1)*rnorm(2)) * c1/cd**3
        dudx(1,3)=dudx(1,3)-(du(1)*rnorm(3)) * c1/cd**3
        dudx(2,1)=dudx(2,1)+(du(2)*rnorm(1)) * c1/cd**3
        dudx(2,2)=dudx(2,2)+(du(2)*rnorm(2)) * c1/cd**3
        dudx(2,3)=dudx(2,3)-(du(2)*rnorm(3)) * c1/cd**3
        dudx(3,1)=dudx(3,1)+(du(3)*rnorm(1)) * c1/cd**3
        dudx(3,2)=dudx(3,2)+(du(3)*rnorm(2)) * c1/cd**3
        dudx(3,3)=dudx(3,3)-(du(3)*rnorm(3)) * c1/cd**3
c
        dudx(1,1)=dudx(1,1)+(du(1)*rnorm(1)) * c1/cd**3
        dudx(1,2)=dudx(1,2)+(du(2)*rnorm(1)) * c1/cd**3
        dudx(1,3)=dudx(1,3)-(du(3)*rnorm(1)) * c1/cd**3
        dudx(2,1)=dudx(2,1)+(du(1)*rnorm(2)) * c1/cd**3
        dudx(2,2)=dudx(2,2)+(du(2)*rnorm(2)) * c1/cd**3
        dudx(2,3)=dudx(2,3)-(du(3)*rnorm(2)) * c1/cd**3
        dudx(3,1)=dudx(3,1)+(du(1)*rnorm(3)) * c1/cd**3
        dudx(3,2)=dudx(3,2)+(du(2)*rnorm(3)) * c1/cd**3
        dudx(3,3)=dudx(3,3)-(du(3)*rnorm(3)) * c1/cd**3
c
        endif
c
c ---------------------------------
c
c       harmonic fix for image source
c
ccc	if (2.ne.3) goto 222
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        ftmp(1)=-dd*dx
        ftmp(2)=-dd*dy
        ftmp(3)=-dd*dz
        fvec(1)=fvec(1)+2*ftmp(1)/cd**3
        fvec(2)=fvec(2)+2*ftmp(2)/cd**3
        fvec(3)=fvec(3)+2*ftmp(3)/cd**3
	if (ifstrain.eq.1) then
          dudx(1,1)=dudx(1,1)+dd*2*(-1/cd**3+ 3*dx*dx/cd**5)
          dudx(1,2)=dudx(1,2)+dd*2*(3*dx*dy/cd**5)
          dudx(1,3)=dudx(1,3)-dd*2*(3*dx*dz/cd**5)
c 
          dudx(2,1)=dudx(2,1)+dd*2*(3*dx*dy/cd**5)
          dudx(2,2)=dudx(2,2)+dd*2*(-1/cd**3+ 3*dy*dy/cd**5)
          dudx(2,3)=dudx(2,3)-dd*2*(3*dy*dz/cd**5)
c
          dudx(3,1)=dudx(3,1)+dd*2*(3*dz*dx/cd**5)
          dudx(3,2)=dudx(3,2)+dd*2*(3*dz*dy/cd**5)
          dudx(3,3)=dudx(3,3)-dd*2*(-1/cd**3+ 3*dz*dz/cd**5)
        endif
222   continue
c ---------------------------------
c
c       ... step 2
c
ccc	if (2.ne.3) goto 333
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
c
        c2=alphaim
ccc        write(6,*) 'rlam,rmu,c2',rlam,rmu,c2
ccc        write(13,*) 'rlam,rmu,c2',rlam,rmu,c2
c
        fvec(1)=fvec(1)+3*c2*dd*dx/cd**5
        fvec(2)=fvec(2)+3*c2*dd*dy/cd**5
        fvec(3)=fvec(3)+3*c2*dd*dz/cd**5
c
        if( ifstrain .eq. 1 ) then
c
        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
c
        dudx(1,1)=dudx(1,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dx/cd**5
        dudx(1,2)=dudx(1,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dx/cd**5
        dudx(1,3)=dudx(1,3)-3*c2*(du(3)*d2+rnorm(3)*d1)*dx/cd**5
        dudx(2,1)=dudx(2,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dy/cd**5
        dudx(2,2)=dudx(2,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dy/cd**5
        dudx(2,3)=dudx(2,3)-3*c2*(du(3)*d2+rnorm(3)*d1)*dy/cd**5
        dudx(3,1)=dudx(3,1)+3*c2*(du(1)*d2+rnorm(1)*d1)*dz/cd**5
        dudx(3,2)=dudx(3,2)+3*c2*(du(2)*d2+rnorm(2)*d1)*dz/cd**5
        dudx(3,3)=dudx(3,3)-3*c2*(du(3)*d2+rnorm(3)*d1)*dz/cd**5
c
        dudx(1,1)=dudx(1,1)+3*c2*dd*(1/cd**5-5*dx*dx/cd**7)
        dudx(1,2)=dudx(1,2)+3*c2*dd*(       -5*dx*dy/cd**7)
        dudx(1,3)=dudx(1,3)-3*c2*dd*(       -5*dx*dz/cd**7)
        dudx(2,1)=dudx(2,1)+3*c2*dd*(       -5*dy*dx/cd**7)
        dudx(2,2)=dudx(2,2)+3*c2*dd*(1/cd**5-5*dy*dy/cd**7)
        dudx(2,3)=dudx(2,3)-3*c2*dd*(       -5*dy*dz/cd**7)
        dudx(3,1)=dudx(3,1)+3*c2*dd*(       -5*dz*dx/cd**7)
        dudx(3,2)=dudx(3,2)+3*c2*dd*(       -5*dz*dy/cd**7)
        dudx(3,3)=dudx(3,3)-3*c2*dd*(1/cd**5-5*dz*dz/cd**7)
      endif
c
333    continue
c
c       ... evaluate strain
c       
      if (ifstrain.eq.1) then
        do i=1,3
        do j=1,3
        strain(i,j)=(dudx(i,j)+dudx(j,i))/2
        enddo
        enddo
      endif
c
ccc      write(6,*)' target at return is ',target
ccc      write(6,*)' source at return is ',source
ccc      write(6,*)' fvec at return is ',fvec
ccc      pause
        return
        end
c
c
c
c
c
