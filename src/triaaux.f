c
c
c
c
c
        subroutine triaarea(vert1,vert2,vert3,area)
        implicit none
        real *8 vert1(3),vert2(3),vert3(3),area
        real *8 x(3),y(3),z(3)
c
c       compute area of an oriented triangle in R^3
c
c       INPUT: 
c   
c       vert1(3),vert2(3),vert3(3)   -  triangle vertices
c
c       OUTPUT: 
c   
c       area          -  triangle area
c
        x(1)=vert2(1)-vert1(1)
        x(2)=vert2(2)-vert1(2)
        x(3)=vert2(3)-vert1(3)
c
        y(1)=vert3(1)-vert1(1)
        y(2)=vert3(2)-vert1(2)
        y(3)=vert3(3)-vert1(3)
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
c
        area=dsqrt(z(1)**2+z(2)**2+z(3)**2)
        area=area*0.5d0
c
        return
        end
c
c
c
c
c
        subroutine trianorm(vert1,vert2,vert3,trinorm)
        implicit none
        real *8 vert1(3),vert2(3),vert3(3),trinorm(3)
        real *8 x(3),y(3),z(3),scale
c
c       construct normal vector of an oriented triangle in R^3
c
c       INPUT: 
c   
c       vert1(3),vert2(3),vert3(3)   -  triangle vertices
c
c       OUTPUT: 
c   
c       trinorm(3)   -  triangle normal
c
c
        x(1)=vert2(1)-vert1(1)
        x(2)=vert2(2)-vert1(2)
        x(3)=vert2(3)-vert1(3)
c
        y(1)=vert3(1)-vert2(1)
        y(2)=vert3(2)-vert2(2)
        y(3)=vert3(3)-vert2(3)
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
c
        scale=dsqrt(z(1)**2+z(2)**2+z(3)**2)
        scale=1/scale
c
        trinorm(1)=z(1)*scale
        trinorm(2)=z(2)*scale
        trinorm(3)=z(3)*scale
c
        return
        end
c
c
c
c
