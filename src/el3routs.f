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
c    $Date: 2012-02-24 09:12:30 -0700 (Fri, 24 Feb 2012) $
c    $Revision: 2736 $
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Elastostatic potential evaluation and decomposition routines     
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine direct(nparts,source,charge,
     1                  point,rpot,ptfrc)
c
c     direct calculation subroutine for electrostatic N-body problem.
c
c     INPUT:
c
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     charge(i) = strength of ith charge
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     rpot =  computed potential
c     ptfrc(3) = computed electrostatic field (-gradient)
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),charge(nparts)
      real *8 point(3),ptfrc(3),rpot
      integer nparts
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c
      do 1000 j = 1,nparts
	 rx = point(1) - source(1,j)
	 ry = point(2) - source(2,j)
	 rz = point(3) - source(3,j)
	 rr = rx*rx + ry*ry + rz*rz
	 rdis = sqrt(rr)
	 rpot = rpot + charge(j)/rdis
	 rmul = charge(j)/(rdis*rr)
	 ptfrc(1) = ptfrc(1) + rmul*rx
	 ptfrc(2) = ptfrc(2) + rmul*ry
	 ptfrc(3) = ptfrc(3) + rmul*rz
1000  continue
      return
      end
c
c
c
c
c
      subroutine directh(nparts,source,charge,
     1                  point,rpot,ptfrc,hessmatr)
c
c     direct calculation subroutine for electrostatic N-body problem,
c     including second deriviatives
c
c     INPUT:
c
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     charge(i) = strength of ith charge
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     rpot =  computed potential
c     ptfrc(3) = computed electrostatic field (-gradient)
c     hessmatr(3,3) =  second derivatives of potential
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),charge(nparts)
      real *8 point(3),ptfrc(3),rpot,hessmatr(3,3)
      integer nparts
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
      do i=1,3
      do j=1,3
         hessmatr(i,j) = 0.0d0
      enddo
      enddo
c
      do 1000 j = 1,nparts
	 rx = point(1) - source(1,j)
	 ry = point(2) - source(2,j)
	 rz = point(3) - source(3,j)
	 rr = rx*rx + ry*ry + rz*rz
	 rdis = sqrt(rr)
	 rpot = rpot + charge(j)/rdis
	 rmul = charge(j)/(rdis*rr)
	 ptfrc(1) = ptfrc(1) + rmul*rx
	 ptfrc(2) = ptfrc(2) + rmul*ry
	 ptfrc(3) = ptfrc(3) + rmul*rz
c
c> r := sqrt(x*x+y*y+z*z);
c                                   2    2    2 1/2
c                            r := (x  + y  + z )
c
c> diff(1/r,x,x);         
c                              2
c                             x                    1
c                    3 ----------------- - -----------------
c                        2    2    2 5/2     2    2    2 3/2
c                      (x  + y  + z )      (x  + y  + z )
c
c> diff(1/r,x,y);
c                                       x y
c                              3 -----------------
c                                  2    2    2 5/2
c                                (x  + y  + z )
c
         rmul2 = charge(j)/(rdis*rr*rr)
         hessmatr(1,1) = hessmatr(1,1) + 3*rx*rx*rmul2 - rmul
         hessmatr(1,2) = hessmatr(1,2) + 3*rx*ry*rmul2
         hessmatr(1,3) = hessmatr(1,3) + 3*rx*rz*rmul2
         hessmatr(2,1) = hessmatr(2,1) + 3*ry*rx*rmul2
         hessmatr(2,2) = hessmatr(2,2) + 3*ry*ry*rmul2 - rmul
         hessmatr(2,3) = hessmatr(2,3) + 3*ry*rz*rmul2
         hessmatr(3,1) = hessmatr(3,1) + 3*rz*rx*rmul2
         hessmatr(3,2) = hessmatr(3,2) + 3*rz*ry*rmul2
         hessmatr(3,3) = hessmatr(3,3) + 3*rz*rz*rmul2 - rmul
1000  continue
      return
      end
c
c
c
c
c
      subroutine directd(nparts,source,charge,dipstr,dipvec,
     1                  point,rpot,ptfrc)
c
c     direct calculation subroutine for electrostatic N-body problem.
c
c     INPUT:
c
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     charge(i) = strength of ith charge
c     dipstr(i) = strength of ith dipole
c     dipvec(3,i) = orientation of ith dipole
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     rpot =  computed potential
c     ptfrc(3) = computed electrostatic field (-gradient)
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),charge(nparts)
      real *8 dipstr(nparts),dipvec(3,nparts)
      real *8 point(3),ptfrc(3),rpot
      integer nparts
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c
      do 1000 j = 1,nparts
         rx = point(1) - source(1,j)
         ry = point(2) - source(2,j)
         rz = point(3) - source(3,j)
         rr = rx*rx + ry*ry + rz*rz
         rdis = sqrt(rr)
         rpot = rpot + charge(j)/rdis
         rmul = charge(j)/(rdis*rr)
         ptfrc(1) = ptfrc(1) + rmul*rx
         ptfrc(2) = ptfrc(2) + rmul*ry
         ptfrc(3) = ptfrc(3) + rmul*rz
         rtmp = dipstr(j)*( dipvec(1,j)*rx
     1               + dipvec(2,j)*ry
     2               + dipvec(3,j)*rz)
         rrr = rdis*rr
         rpot = rpot + rtmp/rrr
         rtmp = 3.0d0*rtmp/(rrr*rr)
         rtmp2 = dipstr(j)/rrr
         ptfrc(1) = ptfrc(1) - dipvec(1,j)*rtmp2
     1              + rtmp*rx
         ptfrc(2) = ptfrc(2) - dipvec(2,j)*rtmp2
     1              + rtmp*ry
         ptfrc(3) = ptfrc(3) - dipvec(3,j)*rtmp2
     1              + rtmp*rz
1000  continue
      return
      end
c
c
c
c
c
      subroutine directdh(nparts,source,charge,dipstr,dipvec,
     1                  point,rpot,ptfrc,hessmatr)
c
c     direct calculation subroutine for electrostatic N-body problem,
c     with second derivatives.
c
c     INPUT:
c
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     charge(i) = strength of ith charge
c     dipstr(i) = strength of ith dipole
c     dipvec(3,i) = orientation of ith dipole
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     rpot =  computed potential
c     ptfrc(3) = computed electrostatic field (-gradient)
c     hessmatr(3,3) =  second derivatives of potential
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),charge(nparts)
      real *8 dipstr(nparts),dipvec(3,nparts)
      real *8 point(3),ptfrc(3),rpot,hessmatr(3,3)
      integer nparts
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
      do i=1,3
      do j=1,3
         hessmatr(i,j) = 0.0d0
      enddo
      enddo
c
      do 1000 j = 1,nparts
         rx = point(1) - source(1,j)
         ry = point(2) - source(2,j)
         rz = point(3) - source(3,j)
         rr = rx*rx + ry*ry + rz*rz
         rdis = sqrt(rr)
         rpot = rpot + charge(j)/rdis
         rmul = charge(j)/(rdis*rr)
         ptfrc(1) = ptfrc(1) + rmul*rx
         ptfrc(2) = ptfrc(2) + rmul*ry
         ptfrc(3) = ptfrc(3) + rmul*rz

c
c> r := sqrt(x*x+y*y+z*z);
c                                   2    2    2 1/2
c                            r := (x  + y  + z )
c
c> diff(1/r,x,x);         
c                              2
c                             x                    1
c                    3 ----------------- - -----------------
c                        2    2    2 5/2     2    2    2 3/2
c                      (x  + y  + z )      (x  + y  + z )
c
c> diff(1/r,x,y);
c                                       x y
c                              3 -----------------
c                                  2    2    2 5/2
c                                (x  + y  + z )
c

         rmul2 = charge(j)/(rdis*rr*rr)
         hessmatr(1,1) = hessmatr(1,1) + 3*rx*rx*rmul2 - rmul
         hessmatr(1,2) = hessmatr(1,2) + 3*rx*ry*rmul2
         hessmatr(1,3) = hessmatr(1,3) + 3*rx*rz*rmul2
         hessmatr(2,1) = hessmatr(2,1) + 3*ry*rx*rmul2
         hessmatr(2,2) = hessmatr(2,2) + 3*ry*ry*rmul2 - rmul
         hessmatr(2,3) = hessmatr(2,3) + 3*ry*rz*rmul2
         hessmatr(3,1) = hessmatr(3,1) + 3*rz*rx*rmul2
         hessmatr(3,2) = hessmatr(3,2) + 3*rz*ry*rmul2
         hessmatr(3,3) = hessmatr(3,3) + 3*rz*rz*rmul2 - rmul

         rtmp = dipstr(j)*( dipvec(1,j)*rx
     1               + dipvec(2,j)*ry
     2               + dipvec(3,j)*rz)
         rrr = rdis*rr
         rpot = rpot + rtmp/rrr
         rtmp = 3.0d0*rtmp/(rrr*rr)
         rtmp2 = dipstr(j)/rrr
         ptfrc(1) = ptfrc(1) - dipvec(1,j)*rtmp2
     1              + rtmp*rx
         ptfrc(2) = ptfrc(2) - dipvec(2,j)*rtmp2
     1              + rtmp*ry
         ptfrc(3) = ptfrc(3) - dipvec(3,j)*rtmp2
     1              + rtmp*rz

c
c> r := sqrt(x*x+y*y+z*z);
c                                   2    2    2 1/2
c                            r := (x  + y  + z )
c
c> diff((a1*x+a2*y+a3*z)/r^3,x,x);
c                                                   2
c            a1 x             (a1 x + a2 y + a3 z) x      a1 x + a2 y + a3 z
c   -6 ----------------- + 15 ----------------------- - 3 ------------------
c        2    2    2 5/2           2    2    2 7/2          2    2    2 5/2
c      (x  + y  + z )            (x  + y  + z )           (x  + y  + z )
c
c
c> diff((a1*x+a2*y+a3*z)/r^3,x,y);
c            a1 y                  a2 x             (a1 x + a2 y + a3 z) x y
c   -3 ----------------- - 3 ----------------- + 15 ------------------------
c        2    2    2 5/2       2    2    2 5/2           2    2    2 7/2
c      (x  + y  + z )        (x  + y  + z )            (x  + y  + z )
c
         rtmp = dipstr(j)*( dipvec(1,j)*rx
     1               + dipvec(2,j)*ry
     2               + dipvec(3,j)*rz)

         hessmatr(1,1)=hessmatr(1,1)+rtmp*3*(5*rx*rx/rdis**7-1/rdis**5)
         hessmatr(1,2)=hessmatr(1,2)+rtmp*3*(5*rx*ry/rdis**7          )
         hessmatr(1,3)=hessmatr(1,3)+rtmp*3*(5*rx*rz/rdis**7          )
         hessmatr(2,1)=hessmatr(2,1)+rtmp*3*(5*ry*rx/rdis**7          )
         hessmatr(2,2)=hessmatr(2,2)+rtmp*3*(5*ry*ry/rdis**7-1/rdis**5)
         hessmatr(2,3)=hessmatr(2,3)+rtmp*3*(5*ry*rz/rdis**7          )
         hessmatr(3,1)=hessmatr(3,1)+rtmp*3*(5*rz*rx/rdis**7          )
         hessmatr(3,2)=hessmatr(3,2)+rtmp*3*(5*rz*ry/rdis**7          )
         hessmatr(3,3)=hessmatr(3,3)+rtmp*3*(5*rz*rz/rdis**7-1/rdis**5)

         hessmatr(1,1)=hessmatr(1,1) 
     $      - 3*dipstr(j)*(dipvec(1,j)*rx+dipvec(1,j)*rx)/rdis**5 
         hessmatr(1,2)=hessmatr(1,2) 
     $      - 3*dipstr(j)*(dipvec(2,j)*rx+dipvec(1,j)*ry)/rdis**5 
         hessmatr(1,3)=hessmatr(1,3) 
     $      - 3*dipstr(j)*(dipvec(3,j)*rx+dipvec(1,j)*rz)/rdis**5 

         hessmatr(2,1)=hessmatr(2,1) 
     $      - 3*dipstr(j)*(dipvec(1,j)*ry+dipvec(2,j)*rx)/rdis**5 
         hessmatr(2,2)=hessmatr(2,2) 
     $      - 3*dipstr(j)*(dipvec(2,j)*ry+dipvec(2,j)*ry)/rdis**5 
         hessmatr(2,3)=hessmatr(2,3) 
     $      - 3*dipstr(j)*(dipvec(3,j)*ry+dipvec(2,j)*rz)/rdis**5 

         hessmatr(3,1)=hessmatr(3,1) 
     $      - 3*dipstr(j)*(dipvec(1,j)*rz+dipvec(3,j)*rx)/rdis**5 
         hessmatr(3,2)=hessmatr(3,2) 
     $      - 3*dipstr(j)*(dipvec(2,j)*rz+dipvec(3,j)*ry)/rdis**5 
         hessmatr(3,3)=hessmatr(3,3) 
     $      - 3*dipstr(j)*(dipvec(3,j)*rz+dipvec(3,j)*rz)/rdis**5 

1000  continue
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine eludirect(rlam,rmu,nparts,source,sigma,
     1                  point,ptfrc)
c
c     direct calculation subroutine for elastostatic 
c     single layer N-body problem.
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     sigma(3,n) = vector strength of nth charge 
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     ptfrc(3) = displacement at the target
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 rvec(3)
      real *8 point(3),ptfrc(3),ptfrc0(3)
      integer nparts
c
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c
      do 1000 j = 1,nparts
	 rx = point(1) - source(1,j)
	 ry = point(2) - source(2,j)
	 rz = point(3) - source(3,j)
	 rr = rx*rx + ry*ry + rz*rz
	 rdis = sqrt(rr)
         rvec(1) = rx
         rvec(2) = ry
         rvec(3) = rz
ccc         call prin2(' rvec is *',rvec,3)
c
         call green3elu(rlam,rmu,rvec,sigma(1,j),ptfrc0)
c
         do i = 1,3
            ptfrc(i) = ptfrc(i) + ptfrc0(i)
         enddo
c
1000  continue
      return
      end
c
c
c
c
c
      subroutine elustdirect(rlam,rmu,nparts,source,sigma,
     1                  point,ptfrc,strain)
c
c     direct calculation subroutine for elastostatic 
c     single layer N-body problem, including strain.
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     sigma(3,n) = vector strength of nth charge
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     ptfrc(3) = displacement at the target
c     strain(3,3) = strain at the target
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 rvec(3)
      real *8 point(3),ptfrc(3),ptfrc0(3)
      real *8 strain(3,3),strain0(3,3)
      integer nparts
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
      do 1000 j = 1,nparts
	 rx = point(1) - source(1,j)
	 ry = point(2) - source(2,j)
	 rz = point(3) - source(3,j)
	 rr = rx*rx + ry*ry + rz*rz
	 rdis = sqrt(rr)
         rvec(1) = rx
         rvec(2) = ry
         rvec(3) = rz
ccc         call prin2(' rvec is *',rvec,3)
c
         call green3elu_strain(rlam,rmu,rvec,sigma(1,j),ptfrc0,strain0)
c
         do i = 1,3
            ptfrc(i) = ptfrc(i) + ptfrc0(i)
         enddo
c
         do i = 1,3
         do k = 1,3
            strain(i,k) = strain(i,k) + strain0(i,k)
         enddo
         enddo
c
1000  continue
      return
      end
C
c
c
c
c
      subroutine elufromlap(rlam,rmu,nparts,source,sigma,
     1                  point,ptfrc)
c
c     direct calculation subroutine for elastostatic 
c     single layer N-body problem, computed from a sequence of 
c     electrostatic interactions.
c
c     INPUT:
c
c     rlam, rmu = Lame parameters
c     nparts = number of sources
c     source(1,i) = x-coordinate of ith charge
c     source(2,i) = y-coordinate of ith charge
c     source(3,i) = z-coordinate of ith charge
c     sigma(3,n) = vector strength of nth charge
c     point(3) = evaluation point
c
c     OUTPUT:
c
c     ptfrc(3) = displacement at the target
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 temp1(3)
      real *8 temp2(3)
      real *8 point(3),ptfrc(3)
      integer nparts
c
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c
      alpha=(rlam+rmu)/(rlam+2*rmu)
c
      do j = 1,3
         do n = 1,nparts
            charge(n) = sigma(j,n)  
         enddo
         call direct(nparts,source,charge,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + rpot *(2-alpha)
         ptfrc(1) = ptfrc(1) + point(j)*temp1(1) *(alpha) 
         ptfrc(2) = ptfrc(2) + point(j)*temp1(2) *(alpha)
         ptfrc(3) = ptfrc(3) + point(j)*temp1(3) *(alpha)
      enddo
ccc      call prin2(' ptfrc is *',ptfrc,3)
      do n = 1,nparts
         charge(n) = sigma(1,n)*source(1,n)  
         charge(n) = charge(n)+sigma(2,n)*source(2,n)  
         charge(n) = charge(n)+sigma(3,n)*source(3,n)  
      enddo
      call direct(nparts,source,charge,
     1                  point,rpot,temp2)
      ptfrc(1) = ptfrc(1) - temp2(1) *(alpha) 
      ptfrc(2) = ptfrc(2) - temp2(2) *(alpha) 
      ptfrc(3) = ptfrc(3) - temp2(3) *(alpha) 
C
      ptfrc(1) = ptfrc(1) / (2*rmu)
      ptfrc(2) = ptfrc(2) / (2*rmu)
      ptfrc(3) = ptfrc(3) / (2*rmu)
C
      return
      end
C
c
c
c
c
      subroutine elufromlap3(rlam,rmu,nparts,source,sigma,
     1                  point,ptfrc)
C
C     direct calculation subroutine for elastostatic 
C     single layer N-body problem, computed from a sequence of 
C     electrostatic interactions.
C
c     Optimized: 4 FMM calls, compatible with eltfromlap3
c
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith charge
C     source(2,i) = y-coordinate of ith charge
C     source(3,i) = z-coordinate of ith charge
C     sigma(3,n) = vector strength of nth charge
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 temp1(3)
      real *8 temp2(3)
      real *8 point(3),ptfrc(3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C
        alpha=(rlam+rmu)/(rlam+2*rmu)
c
      do j = 1,3
         do n = 1,nparts
            charge(n) = sigma(j,n)  / (2*rmu)
         enddo
         call direct(nparts,source,charge,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + rpot *(2-alpha)
         ptfrc(1) = ptfrc(1) + point(j)*temp1(1) *(alpha) 
         ptfrc(2) = ptfrc(2) + point(j)*temp1(2) *(alpha)
         ptfrc(3) = ptfrc(3) + point(j)*temp1(3) *(alpha)
      enddo
ccc      call prin2(' ptfrc is *',ptfrc,3)
      do n = 1,nparts
         charge(n) = sigma(1,n)*source(1,n)  
         charge(n) = charge(n)+sigma(2,n)*source(2,n)  
         charge(n) = charge(n)+sigma(3,n)*source(3,n)  
         charge(n) = charge(n) *(alpha) / (2*rmu)
      enddo
      call direct(nparts,source,charge,
     1                  point,rpot,temp2)
      ptfrc(1) = ptfrc(1) - temp2(1) 
      ptfrc(2) = ptfrc(2) - temp2(2) 
      ptfrc(3) = ptfrc(3) - temp2(3) 
C
      return
      end
C
c
c
c
c
      subroutine elustfromlap(rlam,rmu,nparts,source,sigma,
     1                  point,ptfrc,strain)
c
C     direct calculation subroutine for elastostatic 
c     single layer N-body problem, computed from a sequence of 
c     electrostatic interactions.
C
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith charge
C     source(2,i) = y-coordinate of ith charge
C     source(3,i) = z-coordinate of ith charge
C     sigma(3,n) = vector strength of nth charge
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C     strain(3,3) = strain at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 temp1(3)
      real *8 temp2(3,3)
      real *8 point(3),ptfrc(3),hessmatr(3,3),strain(3,3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
      do i=1,3
      do j=1,3
         hessmatr(i,j) = 0.0d0
      enddo
      enddo
C
        alpha=(rlam+rmu)/(rlam+2*rmu)
c
      do j = 1,3
         do n = 1,nparts
            charge(n) = sigma(j,n)  
         enddo
         call directh(nparts,source,charge,
     1                  point,rpot,temp1,temp2)
         ptfrc(j) = ptfrc(j) + rpot *(2-alpha)
         ptfrc(1) = ptfrc(1) + point(j)*temp1(1) *(alpha) 
         ptfrc(2) = ptfrc(2) + point(j)*temp1(2) *(alpha)
         ptfrc(3) = ptfrc(3) + point(j)*temp1(3) *(alpha)
         hessmatr(j,1) = hessmatr(j,1) - temp1(1) *(2-alpha)
         hessmatr(j,2) = hessmatr(j,2) - temp1(2) *(2-alpha)
         hessmatr(j,3) = hessmatr(j,3) - temp1(3) *(2-alpha)
         hessmatr(1,j) = hessmatr(1,j) + temp1(1) *(alpha)
         hessmatr(2,j) = hessmatr(2,j) + temp1(2) *(alpha)
         hessmatr(3,j) = hessmatr(3,j) + temp1(3) *(alpha)
         hessmatr(1,1) = hessmatr(1,1) - point(j)*temp2(1,1) *(alpha)
         hessmatr(1,2) = hessmatr(1,2) - point(j)*temp2(1,2) *(alpha)
         hessmatr(1,3) = hessmatr(1,3) - point(j)*temp2(1,3) *(alpha)
         hessmatr(2,1) = hessmatr(2,1) - point(j)*temp2(2,1) *(alpha)
         hessmatr(2,2) = hessmatr(2,2) - point(j)*temp2(2,2) *(alpha)
         hessmatr(2,3) = hessmatr(2,3) - point(j)*temp2(2,3) *(alpha)
         hessmatr(3,1) = hessmatr(3,1) - point(j)*temp2(3,1) *(alpha)
         hessmatr(3,2) = hessmatr(3,2) - point(j)*temp2(3,2) *(alpha)
         hessmatr(3,3) = hessmatr(3,3) - point(j)*temp2(3,3) *(alpha)
      enddo
ccc      call prin2(' ptfrc is *',ptfrc,3)
      do n = 1,nparts
         charge(n) = sigma(1,n)*source(1,n)  
         charge(n) = charge(n)+sigma(2,n)*source(2,n)  
         charge(n) = charge(n)+sigma(3,n)*source(3,n)  
      enddo
      call directh(nparts,source,charge,
     1                  point,rpot,temp1,temp2)
      ptfrc(1) = ptfrc(1) - temp1(1) *(alpha) 
      ptfrc(2) = ptfrc(2) - temp1(2) *(alpha) 
      ptfrc(3) = ptfrc(3) - temp1(3) *(alpha) 
      hessmatr(1,1) = hessmatr(1,1) + temp2(1,1) *(alpha)
      hessmatr(1,2) = hessmatr(1,2) + temp2(1,2) *(alpha)
      hessmatr(1,3) = hessmatr(1,3) + temp2(1,3) *(alpha)
      hessmatr(2,1) = hessmatr(2,1) + temp2(2,1) *(alpha)
      hessmatr(2,2) = hessmatr(2,2) + temp2(2,2) *(alpha)
      hessmatr(2,3) = hessmatr(2,3) + temp2(2,3) *(alpha)
      hessmatr(3,1) = hessmatr(3,1) + temp2(3,1) *(alpha)
      hessmatr(3,2) = hessmatr(3,2) + temp2(3,2) *(alpha)
      hessmatr(3,3) = hessmatr(3,3) + temp2(3,3) *(alpha)
C
      ptfrc(1) = ptfrc(1) / (2*rmu)
      ptfrc(2) = ptfrc(2) / (2*rmu)
      ptfrc(3) = ptfrc(3) / (2*rmu)
C
      do i=1,3
      do j=1,3
        strain(i,j)=(hessmatr(i,j)+hessmatr(j,i))/2
        strain(i,j)=strain(i,j)/(2*rmu)
      enddo
      enddo
c
      return
      end
C
C
c
c
C
      subroutine eltdirect(rlam,rmu,nparts,source,sigma,sourcenorm,
     1                  point,ptfrc)
c
c     direct calculation subroutine for elastostatic 
c     double layer N-body problem.
C
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith dipole
C     source(2,i) = y-coordinate of ith dipole
C     source(3,i) = z-coordinate of ith dipole
C     sigma(3,n) = vector strength of nth dipole
C     sourcenorm(3,n) =  orientation vector of nth dipole
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 sourcenorm(3,nparts)
      real *8 rvec(3),tnum
      real *8 point(3),ptfrc(3)
      real *8 ptfrc0(3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C
      do 1000 j = 1,nparts
	 rx = point(1) - source(1,j)
	 ry = point(2) - source(2,j)
	 rz = point(3) - source(3,j)
	 rr = rx*rx + ry*ry + rz*rz
	 rdis = sqrt(rr)
         rvec(1) = rx
         rvec(2) = ry
         rvec(3) = rz
Ccc         call prin2(' rvec is *',rvec,3)
c
         call green3elt(rlam,rmu,rvec,sigma(1,j),sourcenorm(1,j),ptfrc0)
         ptfrc(1) = ptfrc(1) + ptfrc0(1)
         ptfrc(2) = ptfrc(2) + ptfrc0(2)
         ptfrc(3) = ptfrc(3) + ptfrc0(3)
1000  continue
      return
      end
c

C
      subroutine eltstdirect(rlam,rmu,nparts,source,sigma,sourcenorm,
     1                  point,ptfrc,strain)
c
C     direct calculation subroutine for elastostatic 
c     double layer N-body problem, including strain.
C
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith dipole
C     source(2,i) = y-coordinate of ith dipole
C     source(3,i) = z-coordinate of ith dipole
C     sigma(3,n) = vector strength of nth dipole
C     sourcenorm(3,n) =  orientation vector of nth dipole
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C     strain(3,3) = strain at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 sourcenorm(3,nparts)
      real *8 rvec(3),tnum
      real *8 point(3),ptfrc(3),strain(3,3)
      real *8 ptfrc0(3),strain0(3,3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C
      do i=1,3
      do j=1,3
         strain(i,j) = 0.0d0
      enddo
      enddo
c
      do 1000 j = 1,nparts
	 rx = point(1) - source(1,j)
	 ry = point(2) - source(2,j)
	 rz = point(3) - source(3,j)
	 rr = rx*rx + ry*ry + rz*rz
	 rdis = sqrt(rr)
         rvec(1) = rx
         rvec(2) = ry
         rvec(3) = rz
Ccc         call prin2(' rvec is *',rvec,3)
c
         call green3elt_strain
     $      (rlam,rmu,rvec,sigma(1,j),sourcenorm(1,j),ptfrc0,strain0)
         ptfrc(1) = ptfrc(1) + ptfrc0(1)
         ptfrc(2) = ptfrc(2) + ptfrc0(2)
         ptfrc(3) = ptfrc(3) + ptfrc0(3)

         do i = 1,3
         do k = 1,3
            strain(i,k) = strain(i,k) + strain0(i,k)
         enddo
         enddo
C
1000  continue
      return
      end
c

C
      subroutine eltfromlap(rlam,rmu,nparts,source,sigma,sourcenorm,
     1                  point,ptfrc)
c
C     direct calculation subroutine for elastostatic double 
c     layer N-body problem, computed from sequence of electrostatic
c     interactions.
C
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith dipole
C     source(2,i) = y-coordinate of ith dipole
C     source(3,i) = z-coordinate of ith dipole
C     sigma(3,n) = vector strength of nth dipole
C     sourcenorm(3,n) =  orientation vector of nth dipole
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 sourcenorm(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 dipstr(1000000)
      real *8 dipvec(3,1000000)
      real *8 temp1(3)
      real *8 temp2(3)
      real *8 point(3),ptfrc(3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C
      c1 = 3*(rlam+rmu)/(rlam+2*rmu)
      c2 = (rmu)/(rlam+2*rmu)
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0D0
            dipstr(n) = sigma(j,n)
            dipvec(1,n) = sourcenorm(1,n)
            dipvec(2,n) = sourcenorm(2,n)
            dipvec(3,n) = sourcenorm(3,n)
         enddo
         call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + c1*rpot + 6*c2*rpot
         ptfrc(1) = ptfrc(1) + c1*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + c1*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + c1*point(j)*temp1(3)
      enddo
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0D0
            dipstr(n) = sourcenorm(j,n)
            dipvec(1,n) = sigma(1,n)
            dipvec(2,n) = sigma(2,n)
            dipvec(3,n) = sigma(3,n)
         enddo
         call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + c1*rpot + 6*c2*rpot
         ptfrc(1) = ptfrc(1) + c1*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + c1*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + c1*point(j)*temp1(3)
      enddo
Ccc      call prin2(' ptfrc is *',ptfrc,3)
      do n = 1,nparts
         charge(n) = sigma(1,n)*sourcenorm(1,n) + 
     1         sigma(2,n)*sourcenorm(2,n) + 
     2         sigma(3,n)*sourcenorm(3,n) 
         charge(n) = 6*c2*charge(n)
         dipstr(n) = sigma(1,n)*source(1,n) + 
     1         sigma(2,n)*source(2,n) + 
     2         sigma(3,n)*source(3,n) 
         dipstr(n) = c1*dipstr(n)
         dipvec(1,n) = sourcenorm(1,n)
         dipvec(2,n) = sourcenorm(2,n)
         dipvec(3,n) = sourcenorm(3,n)
      enddo
      call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
      ptfrc(1) = ptfrc(1) - temp1(1)
      ptfrc(2) = ptfrc(2) - temp1(2)
      ptfrc(3) = ptfrc(3) - temp1(3)
C
      do n = 1,nparts
         charge(n) = 0.0D0
         dipstr(n) = sourcenorm(1,n)*source(1,n) + 
     1         sourcenorm(2,n)*source(2,n) + 
     2         sourcenorm(3,n)*source(3,n) 
         dipvec(1,n) = sigma(1,n)
         dipvec(2,n) = sigma(2,n)
         dipvec(3,n) = sigma(3,n)
      enddo
      call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
      ptfrc(1) = ptfrc(1) - c1*temp1(1)
      ptfrc(2) = ptfrc(2) - c1*temp1(2)
      ptfrc(3) = ptfrc(3) - c1*temp1(3)
      ptfrc(1) = ptfrc(1)/6
      ptfrc(2) = ptfrc(2)/6
      ptfrc(3) = ptfrc(3)/6
      return
      end
C
c
c
c
c
      subroutine eltfromlap2(rlam,rmu,nparts,source,sigma,sourcenorm,
     1                  point,ptfrc)
c
C     direct calculation subroutine for elastostatic double 
c     layer N-body problem, computed from sequence of electrostatic
c     interactions.
C
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith dipole
C     source(2,i) = y-coordinate of ith dipole
C     source(3,i) = z-coordinate of ith dipole
C     sigma(3,n) = vector strength of nth dipole
C     sourcenorm(3,n) =  orientation vector of nth dipole
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 sourcenorm(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 dipstr(1000000)
      real *8 dipvec(3,1000000)
      real *8 temp1(3)
      real *8 temp2(3)
      real *8 point(3),ptfrc(3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C      
      c1 = 3*(rlam+rmu)/(rlam+2*rmu)
      c2 = (rmu)/(rlam+2*rmu)
c
      alpha = (rlam+rmu)/(rlam+2*rmu)
c
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0D0
            dipstr(n) = sigma(j,n)
            dipvec(1,n) = sourcenorm(1,n)
            dipvec(2,n) = sourcenorm(2,n)
            dipvec(3,n) = sourcenorm(3,n)
         enddo
         call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + 3*(2-alpha)*rpot 
         ptfrc(1) = ptfrc(1) + 3*alpha*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + 3*alpha*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + 3*alpha*point(j)*temp1(3)
      enddo
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0D0
            dipstr(n) = sourcenorm(j,n)
            dipvec(1,n) = sigma(1,n)
            dipvec(2,n) = sigma(2,n)
            dipvec(3,n) = sigma(3,n)
         enddo
         call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + 3*(2-alpha)*rpot
         ptfrc(1) = ptfrc(1) + 3*alpha*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + 3*alpha*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + 3*alpha*point(j)*temp1(3)
      enddo
Ccc      call prin2(' ptfrc is *',ptfrc,3)
      do n = 1,nparts
         charge(n) = sigma(1,n)*sourcenorm(1,n) + 
     1         sigma(2,n)*sourcenorm(2,n) + 
     2         sigma(3,n)*sourcenorm(3,n) 
         charge(n) = 6*c2*charge(n)
         dipstr(n) = sigma(1,n)*source(1,n) + 
     1         sigma(2,n)*source(2,n) + 
     2         sigma(3,n)*source(3,n) 
         dipstr(n) = 3*alpha*dipstr(n)
         dipvec(1,n) = sourcenorm(1,n)
         dipvec(2,n) = sourcenorm(2,n)
         dipvec(3,n) = sourcenorm(3,n)
      enddo
      call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
      ptfrc(1) = ptfrc(1) - temp1(1)
      ptfrc(2) = ptfrc(2) - temp1(2)
      ptfrc(3) = ptfrc(3) - temp1(3)
c
      do n = 1,nparts
         charge(n) = 0.0d0
         dipstr(n) = sourcenorm(1,n)*source(1,n) + 
     1         sourcenorm(2,n)*source(2,n) + 
     2         sourcenorm(3,n)*source(3,n) 
         dipstr(n) = 3*alpha*dipstr(n)
         dipvec(1,n) = sigma(1,n)
         dipvec(2,n) = sigma(2,n)
         dipvec(3,n) = sigma(3,n)
      enddo
      call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
      ptfrc(1) = ptfrc(1) - temp1(1)
      ptfrc(2) = ptfrc(2) - temp1(2)
      ptfrc(3) = ptfrc(3) - temp1(3)
      ptfrc(1) = ptfrc(1)/6
      ptfrc(2) = ptfrc(2)/6
      ptfrc(3) = ptfrc(3)/6
      return
      end
C
c
c
c
c
      subroutine eltfromlap3(rlam,rmu,nparts,source,sigma,sourcenorm,
     1                  point,ptfrc)
c
C     direct calculation subroutine for elastostatic double 
c     layer N-body problem, computed from sequence of electrostatic
c     interactions.
C
c     Optimized: 4 FMM calls, compatible with elufromlap3
c
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith dipole
C     source(2,i) = y-coordinate of ith dipole
C     source(3,i) = z-coordinate of ith dipole
C     sigma(3,n) = vector strength of nth dipole
C     sourcenorm(3,n) =  orientation vector of nth dipole
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 sourcenorm(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 dipstr(1000000)
      real *8 dipvec(3,1000000)
      real *8 temp1(3)
      real *8 temp2(3)
      real *8 point(3),ptfrc(3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C      
      c1 = 3*(rlam+rmu)/(rlam+2*rmu)
      c2 = (rmu)/(rlam+2*rmu)
c
      alpha = (rlam+rmu)/(rlam+2*rmu)
c
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0D0
            dipstr(n) = 1
            dipvec(1,n) = sourcenorm(1,n)*sigma(j,n)
            dipvec(2,n) = sourcenorm(2,n)*sigma(j,n)
            dipvec(3,n) = sourcenorm(3,n)*sigma(j,n)
            dipvec(1,n) = dipvec(1,n)+sigma(1,n)*sourcenorm(j,n)
            dipvec(2,n) = dipvec(2,n)+sigma(2,n)*sourcenorm(j,n)
            dipvec(3,n) = dipvec(3,n)+sigma(3,n)*sourcenorm(j,n)
            dipvec(1,n) = dipvec(1,n)/2
            dipvec(2,n) = dipvec(2,n)/2
            dipvec(3,n) = dipvec(3,n)/2
         enddo
         call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
         ptfrc(j) = ptfrc(j) + (2-alpha)*rpot 
         ptfrc(1) = ptfrc(1) + alpha*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + alpha*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + alpha*point(j)*temp1(3)
      enddo
c
      do n = 1,nparts
         charge(n) = sigma(1,n)*sourcenorm(1,n) + 
     1         sigma(2,n)*sourcenorm(2,n) + 
     2         sigma(3,n)*sourcenorm(3,n) 
         charge(n) = c2*charge(n)
         dipstr(n) = 1 
         rval = (sigma(1,n)*source(1,n) + 
     1         sigma(2,n)*source(2,n) + 
     2         sigma(3,n)*source(3,n))*alpha 
         dipvec(1,n) = sourcenorm(1,n)*rval
         dipvec(2,n) = sourcenorm(2,n)*rval
         dipvec(3,n) = sourcenorm(3,n)*rval
         rval = (sourcenorm(1,n)*source(1,n) + 
     1         sourcenorm(2,n)*source(2,n) + 
     2         sourcenorm(3,n)*source(3,n))*alpha
         dipvec(1,n) = dipvec(1,n)+sigma(1,n)*rval
         dipvec(2,n) = dipvec(2,n)+sigma(2,n)*rval
         dipvec(3,n) = dipvec(3,n)+sigma(3,n)*rval
         dipvec(1,n) = dipvec(1,n)/2
         dipvec(2,n) = dipvec(2,n)/2
         dipvec(3,n) = dipvec(3,n)/2
      enddo
      call directd(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1)
      ptfrc(1) = ptfrc(1) - temp1(1)
      ptfrc(2) = ptfrc(2) - temp1(2)
      ptfrc(3) = ptfrc(3) - temp1(3)
C
      return
      end
C
c
c
c
c
      subroutine eltstfromlap(rlam,rmu,nparts,source,sigma,sourcenorm,
     1                  point,ptfrc,strain)
c
C     direct calculation subroutine for elastostatic double 
C     layer N-body problem, computed from sequence of electrostatic
C     interactions.
C
C     INPUT:
C
C     rlam, rmu = Lame parameters
C     nparts = number of sources
C     source(1,i) = x-coordinate of ith dipole
C     source(2,i) = y-coordinate of ith dipole
C     source(3,i) = z-coordinate of ith dipole
C     sigma(3,n) = vector strength of nth dipole
C     sourcenorm(3,n) =  orientation vector of nth dipole
C     point(3) = evaluation point
C
C     OUTPUT:
C
C     ptfrc(3) = displacement at the target
C     strain(3,3) = strain at the target
C
      implicit real *8 (a-h,o-z)
      real *8 source(3,nparts),sigma(3,nparts)
      real *8 sourcenorm(3,nparts)
      real *8 rvec(3)
      real *8 charge(1000000)
      real *8 dipstr(1000000)
      real *8 dipvec(3,1000000)
      real *8 temp1(3)
      real *8 temp2(3,3),hessmatr(3,3)
      real *8 point(3),ptfrc(3),strain(3,3)
      integer nparts
C
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
C
      do i=1,3
      do j=1,3
         hessmatr(I,J) = 0.0D0
      enddo
      enddo

      c1 = 3*(rlam+rmu)/(rlam+2*rmu)
      c2 = (rmu)/(rlam+2*rmu)
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0d0
            dipstr(n) = sigma(j,n)
            dipvec(1,n) = sourcenorm(1,n)
            dipvec(2,n) = sourcenorm(2,n)
            dipvec(3,n) = sourcenorm(3,n)
         enddo
         call directdh(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1,temp2)
         ptfrc(j) = ptfrc(j) + c1*rpot + 6*c2*rpot
         ptfrc(1) = ptfrc(1) + c1*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + c1*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + c1*point(j)*temp1(3)
         hessmatr(j,1) = hessmatr(j,1) - c1*temp1(1) - 6*c2*temp1(1)
         hessmatr(j,2) = hessmatr(j,2) - c1*temp1(2) - 6*c2*temp1(2)
         hessmatr(j,3) = hessmatr(j,3) - c1*temp1(3) - 6*c2*temp1(3)
         hessmatr(1,j) = hessmatr(1,j) + temp1(1) *c1
         hessmatr(2,j) = hessmatr(2,j) + temp1(2) *c1
         hessmatr(3,j) = hessmatr(3,j) + temp1(3) *c1
         hessmatr(1,1) = hessmatr(1,1) - point(j)*temp2(1,1) *c1
         hessmatr(1,2) = hessmatr(1,2) - point(j)*temp2(1,2) *c1
         hessmatr(1,3) = hessmatr(1,3) - point(j)*temp2(1,3) *c1
         hessmatr(2,1) = hessmatr(2,1) - point(j)*temp2(2,1) *c1
         hessmatr(2,2) = hessmatr(2,2) - point(j)*temp2(2,2) *c1
         hessmatr(2,3) = hessmatr(2,3) - point(j)*temp2(2,3) *c1
         hessmatr(3,1) = hessmatr(3,1) - point(j)*temp2(3,1) *c1
         hessmatr(3,2) = hessmatr(3,2) - point(j)*temp2(3,2) *c1
         hessmatr(3,3) = hessmatr(3,3) - point(j)*temp2(3,3) *c1
      enddo
      do j = 1,3
         do n = 1,nparts
            charge(n) = 0.0D0
            dipstr(n) = sourcenorm(j,n)
            dipvec(1,n) = sigma(1,n)
            dipvec(2,n) = sigma(2,n)
            dipvec(3,n) = sigma(3,n)
         enddo
         call directdh(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1,temp2)
         ptfrc(j) = ptfrc(j) + c1*rpot + 6*c2*rpot
         ptfrc(1) = ptfrc(1) + c1*point(j)*temp1(1)
         ptfrc(2) = ptfrc(2) + c1*point(j)*temp1(2)
         ptfrc(3) = ptfrc(3) + c1*point(j)*temp1(3)
         hessmatr(j,1) = hessmatr(j,1) - c1*temp1(1) - 6*c2*temp1(1)
         hessmatr(j,2) = hessmatr(j,2) - c1*temp1(2) - 6*c2*temp1(2)
         hessmatr(j,3) = hessmatr(j,3) - c1*temp1(3) - 6*c2*temp1(3)
         hessmatr(1,j) = hessmatr(1,j) + temp1(1) *c1
         hessmatr(2,j) = hessmatr(2,j) + temp1(2) *c1
         hessmatr(3,j) = hessmatr(3,j) + temp1(3) *c1
         hessmatr(1,1) = hessmatr(1,1) - point(j)*temp2(1,1) *c1
         hessmatr(1,2) = hessmatr(1,2) - point(j)*temp2(1,2) *c1
         hessmatr(1,3) = hessmatr(1,3) - point(j)*temp2(1,3) *c1
         hessmatr(2,1) = hessmatr(2,1) - point(j)*temp2(2,1) *c1
         hessmatr(2,2) = hessmatr(2,2) - point(j)*temp2(2,2) *c1
         hessmatr(2,3) = hessmatr(2,3) - point(j)*temp2(2,3) *c1
         hessmatr(3,1) = hessmatr(3,1) - point(j)*temp2(3,1) *c1
         hessmatr(3,2) = hessmatr(3,2) - point(j)*temp2(3,2) *c1
         hessmatr(3,3) = hessmatr(3,3) - point(j)*temp2(3,3) *c1
      enddo
Ccc      call prin2(' ptfrc is *',ptfrc,3)
      do n = 1,nparts
         charge(n) = sigma(1,n)*sourcenorm(1,n) + 
     1         sigma(2,n)*sourcenorm(2,n) + 
     2         sigma(3,n)*sourcenorm(3,n) 
         charge(n) = 6*c2*charge(n)
         dipstr(n) = sigma(1,n)*source(1,n) + 
     1         sigma(2,n)*source(2,n) + 
     2         sigma(3,n)*source(3,n) 
         dipstr(n) = c1*dipstr(n)
         dipvec(1,n) = sourcenorm(1,n)
         dipvec(2,n) = sourcenorm(2,n)
         dipvec(3,n) = sourcenorm(3,n)
      enddo
      call directdh(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1,temp2)
      ptfrc(1) = ptfrc(1) - temp1(1)
      ptfrc(2) = ptfrc(2) - temp1(2)
      ptfrc(3) = ptfrc(3) - temp1(3)
C
         hessmatr(1,1) = hessmatr(1,1) + temp2(1,1) 
         hessmatr(1,2) = hessmatr(1,2) + temp2(1,2) 
         hessmatr(1,3) = hessmatr(1,3) + temp2(1,3)
         hessmatr(2,1) = hessmatr(2,1) + temp2(2,1)
         hessmatr(2,2) = hessmatr(2,2) + temp2(2,2)
         hessmatr(2,3) = hessmatr(2,3) + temp2(2,3)
         hessmatr(3,1) = hessmatr(3,1) + temp2(3,1)
         hessmatr(3,2) = hessmatr(3,2) + temp2(3,2)
         hessmatr(3,3) = hessmatr(3,3) + temp2(3,3)
C
      do n = 1,nparts
         charge(n) = 0.0D0
         dipstr(n) = sourcenorm(1,n)*source(1,n) + 
     1         sourcenorm(2,n)*source(2,n) + 
     2         sourcenorm(3,n)*source(3,n)
         dipstr(n) = c1*dipstr(n) 
         dipvec(1,n) = sigma(1,n)
         dipvec(2,n) = sigma(2,n)
         dipvec(3,n) = sigma(3,n)
      enddo
      call directdh(nparts,source,charge, dipstr,dipvec,
     1                  point,rpot,temp1,temp2)
      ptfrc(1) = ptfrc(1) - temp1(1)
      ptfrc(2) = ptfrc(2) - temp1(2)
      ptfrc(3) = ptfrc(3) - temp1(3)
c
         hessmatr(1,1) = hessmatr(1,1) + temp2(1,1) 
         hessmatr(1,2) = hessmatr(1,2) + temp2(1,2) 
         hessmatr(1,3) = hessmatr(1,3) + temp2(1,3)
         hessmatr(2,1) = hessmatr(2,1) + temp2(2,1)
         hessmatr(2,2) = hessmatr(2,2) + temp2(2,2)
         hessmatr(2,3) = hessmatr(2,3) + temp2(2,3)
         hessmatr(3,1) = hessmatr(3,1) + temp2(3,1)
         hessmatr(3,2) = hessmatr(3,2) + temp2(3,2)
         hessmatr(3,3) = hessmatr(3,3) + temp2(3,3)
C
      ptfrc(1) = ptfrc(1)/6
      ptfrc(2) = ptfrc(2)/6
      ptfrc(3) = ptfrc(3)/6
C
      do i=1,3
      do j=1,3
        strain(i,j)=(hessmatr(i,j)+hessmatr(j,i))/2
        strain(i,j)=strain(i,j)/6
      enddo
      enddo
c
      return
      end
C
C
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Auxiliary routines for elastostatic FMM
c
c
        subroutine elfmm3dlap1(npts,j,c1,c2,pot,fld,hess,
     $     point,ifptfrc,ptfrc,ifstrain,hessmatr)
        implicit real *8 (a-h,o-z)
        complex *16 pot(1),fld(3,1),hess(6,1)
        real *8 point(3,1),ptfrc(3,1),hessmatr(3,3,1)
        real *8 rpot,temp1(3),temp2(3,3)
c
        if( ifptfrc .eq. 0 .and. ifstrain .eq. 0 ) return
c
        do k=1,npts
c
        rpot=pot(k)
        temp1(1)=fld(1,k)
        temp1(2)=fld(2,k)
        temp1(3)=fld(3,k)
        if( ifstrain .eq. 1 ) then
        temp2(1,1)=hess(1,k)
        temp2(2,2)=hess(2,k)
        temp2(3,3)=hess(3,k)
        temp2(1,2)=hess(4,k)
        temp2(1,3)=hess(5,k)
        temp2(2,3)=hess(6,k)
        temp2(2,1)=hess(4,k)
        temp2(3,1)=hess(5,k)
        temp2(3,2)=hess(6,k)
        endif
c
        if( ifptfrc .eq. 1 ) then
        ptfrc(j,k)=ptfrc(j,k)+(2-c1)*rpot
        ptfrc(1,k)=ptfrc(1,k)+c1*point(j,k)*temp1(1)
        ptfrc(2,k)=ptfrc(2,k)+c1*point(j,k)*temp1(2)
        ptfrc(3,k)=ptfrc(3,k)+c1*point(j,k)*temp1(3)
        endif
        if( ifstrain .eq. 1 ) then
        hessmatr(j,1,k)=hessmatr(j,1,k)-(2-c1)*temp1(1)
        hessmatr(j,2,k)=hessmatr(j,2,k)-(2-c1)*temp1(2)
        hessmatr(j,3,k)=hessmatr(j,3,k)-(2-c1)*temp1(3)
        hessmatr(1,j,k)=hessmatr(1,j,k)+temp1(1)*c1
        hessmatr(2,j,k)=hessmatr(2,j,k)+temp1(2)*c1
        hessmatr(3,j,k)=hessmatr(3,j,k)+temp1(3)*c1
        hessmatr(1,1,k)=hessmatr(1,1,k)-point(j,k)*temp2(1,1)*c1
        hessmatr(1,2,k)=hessmatr(1,2,k)-point(j,k)*temp2(1,2)*c1
        hessmatr(1,3,k)=hessmatr(1,3,k)-point(j,k)*temp2(1,3)*c1
        hessmatr(2,1,k)=hessmatr(2,1,k)-point(j,k)*temp2(2,1)*c1
        hessmatr(2,2,k)=hessmatr(2,2,k)-point(j,k)*temp2(2,2)*c1
        hessmatr(2,3,k)=hessmatr(2,3,k)-point(j,k)*temp2(2,3)*c1
        hessmatr(3,1,k)=hessmatr(3,1,k)-point(j,k)*temp2(3,1)*c1
        hessmatr(3,2,k)=hessmatr(3,2,k)-point(j,k)*temp2(3,2)*c1
        hessmatr(3,3,k)=hessmatr(3,3,k)-point(j,k)*temp2(3,3)*c1
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
        subroutine elfmm3dlap2(npts,pot,fld,hess,
     $     ifptfrc,ptfrc,ifstrain,hessmatr)
        implicit real *8 (a-h,o-z)
        complex *16 pot(1),fld(3,1),hess(6,1)
        real *8 ptfrc(3,1),hessmatr(3,3,1)
        real *8 rpot,temp1(3),temp2(3,3)
c
        if( ifptfrc .eq. 0 .and. ifstrain .eq. 0 ) return
c
        do k=1,npts
c
        rpot=pot(k)
        temp1(1)=fld(1,k)
        temp1(2)=fld(2,k)
        temp1(3)=fld(3,k)
        if( ifstrain .eq. 1 ) then
        temp2(1,1)=hess(1,k)
        temp2(2,2)=hess(2,k)
        temp2(3,3)=hess(3,k)
        temp2(1,2)=hess(4,k)
        temp2(1,3)=hess(5,k)
        temp2(2,3)=hess(6,k)
        temp2(2,1)=hess(4,k)
        temp2(3,1)=hess(5,k)
        temp2(3,2)=hess(6,k)
        endif
c
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,k)=ptfrc(1,k)-temp1(1)
        ptfrc(2,k)=ptfrc(2,k)-temp1(2)
        ptfrc(3,k)=ptfrc(3,k)-temp1(3)
        endif
c
        if( ifstrain .eq. 1 ) then
        hessmatr(1,1,k)=hessmatr(1,1,k)+temp2(1,1) 
        hessmatr(1,2,k)=hessmatr(1,2,k)+temp2(1,2) 
        hessmatr(1,3,k)=hessmatr(1,3,k)+temp2(1,3)
        hessmatr(2,1,k)=hessmatr(2,1,k)+temp2(2,1)
        hessmatr(2,2,k)=hessmatr(2,2,k)+temp2(2,2)
        hessmatr(2,3,k)=hessmatr(2,3,k)+temp2(2,3)
        hessmatr(3,1,k)=hessmatr(3,1,k)+temp2(3,1)
        hessmatr(3,2,k)=hessmatr(3,2,k)+temp2(3,2)
        hessmatr(3,3,k)=hessmatr(3,3,k)+temp2(3,3)
        endif
c
        enddo
c
        return
        end
c
c
c
