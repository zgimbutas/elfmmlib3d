c%function [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
c% CalcTriDisps.m
c%
c% Calculates displacements due to slip on a triangular dislocation in an
c% elastic half space utilizing the Comninou and Dunders (1975) expressions
c% for the displacements due to an angular dislocation in an elastic half
c% space.
c%
c% Arguments
c%  sx : x-coordinates of observation points
c%  sy : y-coordinates of observation points
c%  sz : z-coordinates of observation points
c%  x  : x-coordinates of triangle vertices.
c%  y  : y-coordinates of triangle vertices.
c%  z  : z-coordinates of triangle vertices.
c%  pr : Poisson's ratio
c%  ss : strike slip displacement
c%  ts : tensile slip displacement
c%  ds : dip slip displacement
c%
c% Returns
c%  U  : structure containing the displacements (U.x, U.y, U.z)
c%
c% This paper should and related code should be cited as:
c% Brendan J. Meade, Algorithms for the calculation of exact 
c% displacements, strains, and stresses for Triangular Dislocation 
c% Elements in a uniform elastic half space, Computers & 
c% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
c%
c% Use at your own risk and please let me know of any bugs/errors!
c%
c% Copyright (c) 2006 Brendan Meade
c% 
c% Permission is hereby granted, free of charge, to any person obtaining a
c% copy of this software and associated documentation files (the
c% "Software"), to deal in the Software without restriction, including
c% without limitation the rights to use, copy, modify, merge, publish,
c% distribute, sublicense, and/or sell copies of the Software, and to permit
c% persons to whom the Software is furnished to do so, subject to the
c% following conditions:
c%
c% The above copyright notice and this permission notice shall be included
c% in all copies or substantial portions of the Software.
c% 
c% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
c% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
c% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
c% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
c% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
c% USE OR OTHER DEALINGS IN THE SOFTWARE.


c%function [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds)
c% CalcTriStrains.m
c%
c% Calculates strains due to slip on a triangular dislocation in an
c% elastic half space utilizing the symbolically differentiated
c% displacement gradient tensor derived from the expressions for
c% the displacements due to an angular dislocation in an elastic half
c% space (Comninou and Dunders, 1975).
c%
c% Arguments
c%  sx : x-coordinates of observation points
c%  sy : y-coordinates of observation points
c%  sz : z-coordinates of observation points
c%  x  : x-coordinates of triangle vertices.
c%  y  : y-coordinates of triangle vertices.
c%  z  : z-coordinates of triangle vertices.
c%  pr : Poisson's ratio
c%  ss : strike slip displacement
c%  ts : tensile slip displacement
c%  ds : dip slip displacement
c%
c% Returns
c%  S  : structure containing the strains (S.xx, S.yy, S.zz, S.xy, S.xz, S.yz)
c%
c% This paper should and related code should be cited as:
c% Brendan J. Meade, Algorithms for the calculation of exact 
c% displacements, strains, and stresses for Triangular Dislocation 
c% Elements in a uniform elastic half space, Computers & 
c% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
c%
c% Use at your own risk and please let me know of any bugs/errors.
c%
c% Copyright (c) 2006 Brendan Meade
c% 
c% Permission is hereby granted, free of charge, to any person obtaining a
c% copy of this software and associated documentation files (the
c% "Software"), to deal in the Software without restriction, including
c% without limitation the rights to use, copy, modify, merge, publish,
c% distribute, sublicense, and/or sell copies of the Software, and to permit
c% persons to whom the Software is furnished to do so, subject to the
c% following conditions:
c% 
c% The above copyright notice and this permission notice shall be included
c% in all copies or substantial portions of the Software.
c% 
c% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
c% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
c% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
c% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
c% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
c% USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c
c
c

        real *8 function cot(x)
        implicit real *8 (a-h,o-z)
        cot = 1/tan(x)
        return
        end
c
c
c
c
c
        subroutine  adv(y1, y2, y3, a, beta, nu, B1, B2, B3, 
     $     v1, v2, v3)
        implicit real *8 (a-h,o-z)

c
c These are the displacements in a uniform elastic half space due to slip
c on an angular dislocation (Comninou and Dunders, 1975).  Some of the
c equations for the B2 and B3 cases have been corrected following Thomas
c 1993.  The equations are coded in way such that they roughly correspond
c to each line in original text.  Exceptions have been made where it made 
c more sense because of grouping symbols.
c
        done=1
        pi=4*atan(done)

        sinbeta = sin(beta)
        cosbeta = cos(beta)
        cotbeta = cot(beta)
        z1 = y1*cosbeta - y3*sinbeta
        z3 = y1*sinbeta + y3*cosbeta
        R2 = y1*y1 + y2*y2 + y3*y3
        R = sqrt(R2)
        y3bar = y3 + 2*a
        z1bar = y1*cosbeta + y3bar*sinbeta
        z3bar = -y1*sinbeta + y3bar*cosbeta
        R2bar = y1*y1 + y2*y2 + y3bar*y3bar
        Rbar = sqrt(R2bar)
        F = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2*R *sinbeta, y1*z1
     $     +(y2*y2)*cosbeta)
        Fbar = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2 *Rbar
     $     *sinbeta, y1*z1bar+(y2*y2)*cosbeta)

c
c Case I: Burgers vector (B1,0,0)
c

        v1InfB1 = 2*(1-nu)*(F+Fbar) - y1*y2*(1/(R*(R-y3)) + 1 /(Rbar
     $     *(Rbar+y3bar))) - y2*cosbeta*((R*sinbeta-y1)/(R*(R-z3)) +
     $     (Rbar*sinbeta-y1)/(Rbar*(Rbar+z3bar)))
        v2InfB1 = (1-2*nu)*(log(R-y3)+log(Rbar+y3bar) - cosbeta*(log(R
     $     -z3)+log(Rbar+z3bar))) - y2*y2*(1/(R*(R-y3))+1 /(Rbar*(Rbar
     $     +y3bar)) - cosbeta*(1/(R*(R-z3))+1/(Rbar*(Rbar +z3bar))))
        v3InfB1 = y2 * (1/R - 1/Rbar - cosbeta*((R*cosbeta-y3) /(R*(R
     $     -z3)) - (Rbar*cosbeta+y3bar)/(Rbar*(Rbar+z3bar))))
        v1InfB1 = v1InfB1 / (8*pi*(1-nu))
        v2InfB1 = v2InfB1 / (8*pi*(1-nu))
        v3InfB1 = v3InfB1 / (8*pi*(1-nu))

        v1CB1 = -2*(1-nu)*(1-2*nu)*Fbar*(cotbeta*cotbeta) + (1-2*nu)*y2
     $     /(Rbar+y3bar) * ((1-2*nu-a/Rbar)*cotbeta - y1 /(Rbar+y3bar)
     $     *(nu+a/Rbar)) + (1-2*nu)*y2*cosbeta*cotbeta /(Rbar+z3bar)
     $     *(cosbeta+a/Rbar) + a*y2*(y3bar-a)*cotbeta/(Rbar *Rbar*Rbar)
     $     + y2*(y3bar-a)/(Rbar*(Rbar+y3bar))*(-(1-2*nu) *cotbeta + y1
     $     /(Rbar+y3bar) * (2*nu+a/Rbar) + a*y1/(Rbar *Rbar)) + y2
     $     *(y3bar-a)/(Rbar*(Rbar+z3bar))*(cosbeta/(Rbar +z3bar)*((Rbar
     $     *cosbeta+y3bar) * ((1-2*nu)*cosbeta-a/Rbar) *cotbeta + 2*(1
     $     -nu)*(Rbar*sinbeta-y1)*cosbeta) - a*y3bar *cosbeta*cotbeta
     $     /(Rbar*Rbar))
        v2CB1 = (1-2*nu)*((2*(1-nu)*(cotbeta*cotbeta)-nu) *log(Rbar
     $     +y3bar) -(2*(1-nu)*(cotbeta*cotbeta)+1-2*nu)*cosbeta
     $     *log(Rbar+z3bar)) - (1-2*nu)/(Rbar+y3bar)*(y1*cotbeta*(1-2*nu
     $     -a/Rbar) + nu*y3bar - a + (y2*y2)/(Rbar+y3bar)*(nu+a/Rbar)) -
     $     (1-2*nu)*z1bar*cotbeta/(Rbar+z3bar)*(cosbeta+a/Rbar) - a*y1
     $     *(y3bar-a)*cotbeta/(Rbar*Rbar*Rbar) + (y3bar-a)/(Rbar+y3bar)
     $     *(-2*nu + 1/Rbar*((1-2*nu)*y1*cotbeta-a) + (y2*y2)/(Rbar
     $     *(Rbar+y3bar))*(2*nu+a/Rbar)+a*(y2*y2)/(Rbar*Rbar*Rbar)) +
     $     (y3bar-a)/(Rbar+z3bar)*((cosbeta*cosbeta) - 1/Rbar*((1-2*nu)
     $     *z1bar*cotbeta+a*cosbeta) + a*y3bar*z1bar*cotbeta/(Rbar*Rbar
     $     *Rbar) - 1/(Rbar*(Rbar+z3bar)) * ((y2*y2)*(cosbeta*cosbeta) -
     $     a*z1bar*cotbeta/Rbar*(Rbar*cosbeta+y3bar)))

        v3CB1 = 2*(1-nu)*(((1-2*nu)*Fbar*cotbeta) + (y2 /(Rbar+y3bar)*(2
     $     *nu+a/Rbar)) - (y2*cosbeta/(Rbar+z3bar) *(cosbeta+a/Rbar))) +
     $     y2*(y3bar-a)/Rbar*(2*nu/(Rbar+y3bar)+a /(Rbar*Rbar)) + y2
     $     *(y3bar-a)*cosbeta/(Rbar*(Rbar+z3bar))*(1-2 *nu-(Rbar*cosbeta
     $     +y3bar)/(Rbar+z3bar)*(cosbeta + a/Rbar) - a *y3bar/(Rbar
     $     *Rbar))

        v1CB1 = v1CB1 / (4*pi*(1-nu))
        v2CB1 = v2CB1 / (4*pi*(1-nu))
        v3CB1 = v3CB1 / (4*pi*(1-nu))
        
        v1B1 = v1InfB1 + v1CB1
        v2B1 = v2InfB1 + v2CB1
        v3B1 = v3InfB1 + v3CB1

c
c Case II: Burgers vector (0,B2,0)
c
        v1InfB2 = -(1-2*nu)*(log(R-y3) + log(Rbar+y3bar) -cosbeta*(log(R
     $     -z3)+log(Rbar+z3bar))) + y1*y1*(1/(R*(R-y3))+1 /(Rbar*(Rbar
     $     +y3bar))) + z1*(R*sinbeta-y1)/(R*(R-z3)) + z1bar *(Rbar
     $     *sinbeta -y1)/(Rbar*(Rbar+z3bar))
        v2InfB2 = 2*(1-nu)*(F+Fbar) + y1*y2*(1/(R*(R-y3))+1 /(Rbar*(Rbar
     $     +y3bar))) - y2*(z1/(R*(R-z3))+z1bar/(Rbar*(Rbar +z3bar)))
        v3InfB2 = -(1-2*nu)*sinbeta*(log(R-z3)-log(Rbar +z3bar)) - y1*(1
     $     /R-1/Rbar) + z1*(R*cosbeta-y3)/(R*(R-z3)) - z1bar*(Rbar
     $     *cosbeta+y3bar)/(Rbar*(Rbar+z3bar))
        v1InfB2 = v1InfB2 / (8*pi*(1-nu))
        v2InfB2 = v2InfB2 / (8*pi*(1-nu))
        v3InfB2 = v3InfB2 / (8*pi*(1-nu))

        v1CB2 = (1-2*nu)*((2*(1-nu)*(cotbeta*cotbeta)+nu) *log(Rbar
     $     +y3bar) - (2*(1-nu)*(cotbeta*cotbeta)+1)*cosbeta *log(Rbar
     $     +z3bar)) +(1-2*nu)/(Rbar+y3bar)* (-(1-2*nu)*y1 *cotbeta+nu
     $     *y3bar-a+a*y1*cotbeta/Rbar + (y1*y1)/(Rbar+y3bar) *(nu+a
     $     /Rbar)) - (1-2*nu)*cotbeta/(Rbar+z3bar)*(z1bar*cosbeta - a
     $     *(Rbar*sinbeta-y1)/(Rbar*cosbeta)) - a*y1*(y3bar-a) *cotbeta
     $     /(Rbar*Rbar*Rbar) + (y3bar-a)/(Rbar+y3bar)*(2*nu + 1 /Rbar
     $     *((1-2*nu)*y1*cotbeta+a) - (y1*y1)/(Rbar*(Rbar+y3bar)) *(2*nu
     $     +a/Rbar) - a*(y1*y1)/(Rbar*Rbar*Rbar)) + (y3bar-a) *cotbeta
     $     /(Rbar+z3bar)*(-cosbeta*sinbeta+a*y1*y3bar/(Rbar*Rbar *Rbar
     $     *cosbeta) + (Rbar*sinbeta-y1)/Rbar*(2*(1-nu)*cosbeta - (Rbar
     $     *cosbeta+y3bar)/(Rbar+z3bar)*(1+a/(Rbar*cosbeta))))
        v2CB2 = 2*(1-nu)*(1-2*nu)*Fbar*cotbeta*cotbeta + (1 -2*nu)*y2
     $     /(Rbar+y3bar)*(-(1-2*nu-a/Rbar)*cotbeta + y1/(Rbar +y3bar)
     $     *(nu+a/Rbar)) - (1-2*nu)*y2*cotbeta/(Rbar+z3bar)*(1+a /(Rbar
     $     *cosbeta)) - a*y2*(y3bar-a)*cotbeta/(Rbar*Rbar*Rbar) + y2
     $     *(y3bar-a)/(Rbar*(Rbar+y3bar))*((1-2*nu)*cotbeta - 2*nu*y1
     $     /(Rbar+y3bar) - a*y1/Rbar*(1/Rbar+1/(Rbar+y3bar))) + y2
     $     *(y3bar-a)*cotbeta/(Rbar*(Rbar+z3bar))*(-2*(1-nu)*cosbeta +
     $     (Rbar*cosbeta+y3bar)/(Rbar+z3bar)*(1+a/(Rbar*cosbeta)) + a
     $     *y3bar/((Rbar*Rbar)*cosbeta))
        v3CB2 = -2*(1-nu)*(1-2*nu)*cotbeta * (log(Rbar +y3bar)-cosbeta
     $     *log(Rbar+z3bar)) - 2*(1-nu)*y1/(Rbar+y3bar) *(2*nu+a/Rbar) +
     $     2*(1-nu)*z1bar/(Rbar+z3bar)*(cosbeta+a/Rbar) + (y3bar-a)/Rbar
     $     *((1-2*nu)*cotbeta-2*nu*y1/(Rbar+y3bar)-a*y1 /(Rbar*Rbar)) -
     $     (y3bar-a)/(Rbar+z3bar)*(cosbeta*sinbeta + (Rbar*cosbeta
     $     +y3bar)*cotbeta/Rbar*(2*(1-nu)*cosbeta - (Rbar *cosbeta
     $     +y3bar)/(Rbar+z3bar)) + a/Rbar*(sinbeta - y3bar*z1bar /(Rbar
     $     *Rbar) - z1bar*(Rbar*cosbeta+y3bar)/(Rbar*(Rbar +z3bar))))
        v1CB2 = v1CB2 / (4*pi*(1-nu))
        v2CB2 = v2CB2 / (4*pi*(1-nu))
        v3CB2 = v3CB2 / (4*pi*(1-nu))

        v1B2 = v1InfB2 + v1CB2
        v2B2 = v2InfB2 + v2CB2
        v3B2 = v3InfB2 + v3CB2

c
c Case III: Burgers vector (0,0,B3)
c
        v1InfB3 = y2*sinbeta*((R*sinbeta-y1)/(R*(R-z3))+(Rbar*sinbeta
     $     -y1)/(Rbar*(Rbar+z3bar)))
        v2InfB3 = (1-2*nu)*sinbeta*(log(R-z3)+log(Rbar+z3bar)) - (y2*y2)
     $     *sinbeta*(1/(R*(R-z3))+1/(Rbar*(Rbar+z3bar)))
        v3InfB3 = 2*(1-nu)*(F-Fbar) + y2*sinbeta*((R*cosbeta-y3)/(R*(R
     $     -z3))-(Rbar*cosbeta+y3bar)/(Rbar*(Rbar+z3bar)))

        v1InfB3 = v1InfB3 / (8*pi*(1-nu))
        v2InfB3 = v2InfB3 / (8*pi*(1-nu))
        v3InfB3 = v3InfB3 / (8*pi*(1-nu))

        v1CB = (1-2*nu)*(y2/(Rbar+y3bar)*(1+a/Rbar) - y2*cosbeta /(Rbar
     $     +z3bar)*(cosbeta+a/Rbar)) -y2*(y3bar-a)/Rbar*(a /(Rbar*Rbar)
     $     + 1 /(Rbar+y3bar)) + y2*(y3bar-a)*cosbeta /(Rbar*(Rbar
     $     +z3bar)) *((Rbar*cosbeta+y3bar)/(Rbar+z3bar) *(cosbeta+a
     $     /Rbar) + a*y3bar /(Rbar*Rbar))
        v2CB3 = (1-2*nu)*(-sinbeta*log(Rbar+z3bar) - y1 /(Rbar+y3bar)*(1
     $     +a/Rbar) + z1bar/(Rbar+z3bar)*(cosbeta+a /Rbar)) + y1*(y3bar
     $     -a)/Rbar*(a/(Rbar*Rbar) + 1/(Rbar +y3bar)) - (y3bar-a)/(Rbar
     $     +z3bar)*(sinbeta*(cosbeta-a /Rbar) + z1bar /Rbar*(1+a*y3bar
     $     /(Rbar*Rbar)) - 1/(Rbar *(Rbar+z3bar))*((y2*y2)*cosbeta
     $     *sinbeta - a*z1bar/Rbar *(Rbar*cosbeta+y3bar)))
        v3CB3 = 2*(1-nu)*Fbar + 2*(1-nu)*(y2*sinbeta/(Rbar+z3bar)
     $     *(cosbeta + a/Rbar)) + y2*(y3bar-a)*sinbeta/(Rbar*(Rbar
     $     +z3bar))*(1 + (Rbar*cosbeta+y3bar)/(Rbar+z3bar)*(cosbeta +a
     $     /Rbar) + a*y3bar/(Rbar*Rbar))

        v1CB3 = v1CB3 / (4*pi*(1-nu))
        v2CB3 = v2CB3 / (4*pi*(1-nu))
        v3CB3 = v3CB3 / (4*pi*(1-nu))

        v1B3 = v1InfB3 + v1CB3
        v2B3 = v2InfB3 + v2CB3
        v3B3 = v3InfB3 + v3CB3

c
c Sum the for each slip component
c
        v1 = B1*v1B1 + B2*v1B2 + B3*v1B3
        v2 = B1*v2B1 + B2*v2B2 + B3*v2B3
        v3 = B1*v3B1 + B2*v3B2 + B3*v3B3

        return
        end
c
c
c
c
c
        subroutine advs(y1, y2, y3, a, b, nu, B1, B2, B3,
     $     e11, e22, e33, e12, e13, e23)
        implicit real *8 (a-h,o-z)
c        
c These are the strains in a uniform elastic half space due to slip
c on an angular dislocation.  They were calculated by symbolically
c differentiating the expressions for the displacements (Comninou and
c Dunders, 1975, with typos noted by Thomas 1993) then combining the
c elements of the displacement gradient tensor to form the strain tensor.
c
        done=1
        pi=4*atan(done)
 
        e11 = B1*(0.125d0*((2-2*nu)*(2*y2/y1**2/(1+y2**2/y1**2)-y2/(y1
     $     *cos(b)-y3*sin(b))**2*cos(b)/(1+y2**2/(y1*cos(b)-y3*sin(b))
     $     **2)+(y2/(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)
     $     -y3*sin(b))+y2**2*cos(b))*y1-y2*(y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2*(2*y1
     $     *cos(b)-y3*sin(b)))/(1+y2**2*(y1**2+y2**2+y3**2)*sin(b)**2
     $     /(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2)-y2/(y1*cos(b)
     $     +(y3+2*a)*sin(b))**2*cos(b)/(1+y2**2/(y1*cos(b)+(y3+2*a)
     $     *sin(b))**2)+(y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     /(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))*y1-y2*(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)+(y3+2*a)*sin(b)))/(1
     $     +y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2))-y2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a))-y1*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y3)*y1-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y3)**2*y1-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*y1-1
     $     /(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*y1)-y2*cos(b)*((1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*sin(b)*y1-1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-((y1**2+y2**2+y3
     $     **2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))*y1-((y1**2+y2**2
     $     +y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2
     $     +y2**2+y3**2)**(0.5d0)*y1-sin(b))+(1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*sin(b)*y1-1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1-((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))))/pi/(1-nu)+1/4*((-2+2*nu)*(1-2*nu)*(y2/y1**2/(1+y2
     $     **2/y1**2)-y2/(y1*cos(b)+(y3+2*a)*sin(b))**2*cos(b)/(1+y2**2
     $     /(y1*cos(b)+(y3+2*a)*sin(b))**2)+(y2/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))*y1-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1
     $     *(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)
     $     +(y3+2*a)*sin(b)))/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     **2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)
     $     **2-(1-2*nu)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*((1-2*nu-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)-y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1+(1-2*nu)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)*(a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1*cot(b)-1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))+y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y1**2/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0))-(1-2*nu)*y2*cos(b)*cot(b)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0))*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1-sin(b))-(1-2*nu)*y2*cos(b)*cot(b)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)*y1-3*a*y2*(y3+a)*cot(b)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(5/2)*y1-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*((
     $     -1+2*nu)*cot(b)+y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*y1/(y1**2
     $     +y2**2+(y3+2*a)**2))*y1-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2*nu)
     $     *cot(b)+y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu
     $     +a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*y1/(y1**2+y2**2+(y3
     $     +2*a)**2))*y1+y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))-y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)+a/(y1
     $     **2+y2**2+(y3+2*a)**2)-2*a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **2)-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     *((1-2*nu)*cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     *cot(b)+(2-2*nu)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     -y1)*cos(b))-a*(y3+2*a)*cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)
     $     **2))*y1-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)*((1-2*nu)*cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     *cot(b)+(2-2*nu)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     -y1)*cos(b))-a*(y3+2*a)*cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)
     $     **2))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+y2*(y3
     $     +a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(-cos(b)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)
     $     *cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)
     $     *((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+cos(b)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*y1*((1-2*nu)
     $     *cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*y1*cot(b)+(2-2*nu)*(1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*sin(b)*y1-1)*cos(b))+2*a*(y3+2*a)*cos(b)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**2*y1))/pi/(1-nu))+B2
     $     *(0.125d0*((-1+2*nu)*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     -cos(b)*((1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))+(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))))+2*y1*(1/(y1**2+y2**2
     $     +y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a))+y1**2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y3)*y1-1/(y1**2+y2**2+y3**2)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y3)**2*y1-1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *y1-1/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*y1)+cos(b)*((y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)-y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))+(y1*cos(b)-y3*sin(b))*(1/(y1
     $     **2+y2**2+y3**2)**(0.5d0)*sin(b)*y1-1)/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     -(y1*cos(b)-y3*sin(b))*((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)
     $     -y1)/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))*y1-(y1*cos(b)-y3*sin(b))*((y1
     $     **2+y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     **2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))+cos(b)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))+(y1*cos(b)+(y3+2*a)*sin(b))*(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*sin(b)*y1-1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))-(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*y1-(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b)))/pi/(1-nu)+1/4*((1-2*nu)*(((2-2*nu)*cot(b)**2+nu)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)-((2-2*nu)*cot(b)**2+1)*cos(b)*(1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-(1-2*nu)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2*nu)*y1*cot(b)
     $     +nu*(y3+2*a)-a+a*y1*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1+(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)*((-1+2*nu)*cot(b)+a*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-a*y1**2*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)+2
     $     *y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))-y1**3/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)**2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1**3/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0))+(1-2*nu)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*((y1*cos(b)+(y3+2*a)
     $     *sin(b))*cos(b)-a*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     -y1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))*(1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))-(1-2*nu)*cot(b)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(cos(b)**2-a*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)*y1
     $     -1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b)+a*((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/cos(b)*y1)-a*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)+3*a*y1**2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(5/2)-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)**2*(2*nu+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)
     $     *y1*cot(b)+a)-y1**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))-a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1+(y3+a)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)*((1-2*nu)*y1*cot(b)+a)*y1+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(1-2*nu)*cot(b)-2*y1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+y1**3/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     +y1**3/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+y1**3/(y1**2+y2**2+(y3+2*a)**2)**2/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*a-2*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*y1+3*a*y1**3/(y1**2+y2**2+(y3+2*a)**2)**(5/2))
     $     -(y3+a)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*(-cos(b)*sin(b)+a*y1*(y3+2*a)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)/cos(b)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((2-2
     $     *nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))))*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+(y3+a)*cot(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /cos(b)-3*a*y1**2*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)
     $     /cos(b)+(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)*y1-1)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((2-2*nu)*cos(b)-((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/cos(b)))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((2-2
     $     *nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b)))*y1
     $     +((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))
     $     +((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1
     $     +a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))*(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1-sin(b))+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/cos(b)*y1)))/pi/(1-nu))+B3*(0.125d0*y2*sin(b)*((1
     $     /(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)*y1-1)/(y1**2+y2**2+y3
     $     **2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3
     $     *cos(b))-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2
     $     **2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)
     $     -y3*cos(b))*y1-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1
     $     -sin(b))+(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)*y1-1)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*y1-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1-sin(b)))/pi/(1-nu)+1/4*((1-2*nu)*(
     $     -y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(1+a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1-y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1+y2*cos(b)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+y2*cos(b)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)*y1)+y2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a))*y1-y2*(y3+a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*(-2*a/(y1**2+y2**2+(y3+2*a)**2)**2*y1-1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1)-y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))*y1-y2*(y3+a)
     $     *cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2
     $     *a)**2))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+y2
     $     *(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*y1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1-2*a*(y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**2*y1))/pi/(1-nu))

        e22 = B1*(0.125d0*((1-2*nu)*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y2
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     -cos(b)*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y2/((y1**2+y2**2+y3
     $     **2)**(0.5d0)-y1*sin(b)-y3*cos(b))+1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))))-2*y2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)-cos(b)*(1/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3
     $     **2)**(0.5d0)-y1*sin(b)-y3*cos(b))+1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))))-y2**2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y3)*y2-1/(y1**2+y2**2+y3**2)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y3)**2*y2-1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *y2-1/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*y2-cos(b)*(-1/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y2-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))**2*y2-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*y2-1/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*y2)))/pi/(1
     $     -nu)+1/4*((1-2*nu)*(((2-2*nu)*cot(b)**2-nu)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)-((2-2*nu)*cot(b)**2+1-2*nu)*cos(b)/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b)))+(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(y1*cot(b)*(1-2*nu-a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))+nu*(y3+2*a)-a+y2**2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-(1-2*nu)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(a*y1*cot(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)*y2+2*y2/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))-y2**3/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y2**3/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0))+(1-2
     $     *nu)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2+(1-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2+3*a*y2*(y3+a)*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(5/2)*y1-(y3+a)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)**2*(-2*nu+1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*((1-2*nu)*y1*cot(b)-a)+y2**2/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*y2**2/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((1-2*nu)*y1
     $     *cot(b)-a)*y2+2*y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))-y2**3/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))-y2**3/(y1**2+y2**2+(y3+2*a)**2)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(2*nu+a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0))-y2**3/(y1**2+y2**2+(y3+2*a)
     $     **2)**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a+2*a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)*y2-3*a*y2**3/(y1**2+y2**2+(y3
     $     +2*a)**2)**(5/2))-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)**2-1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*((1-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))
     $     *cot(b)+a*cos(b))+a*(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)**2-a*(y1*cos(b)+(y3+2
     $     *a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y2+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*((1-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)+a
     $     *cos(b))*y2-3*a*(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(5/2)*y2+1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(y2**2*cos(b)**2-a*(y1*cos(b)+(y3+2*a)
     $     *sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*y2+1/(y1**2+y2**2
     $     +(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*(y2**2*cos(b)**2-a*(y1*cos(b)+(y3+2*a)
     $     *sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*y2-1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(2*y2*cos(b)**2+a*(y1*cos(b)+(y3+2
     $     *a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*y2-a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)*cos(b)
     $     *y2)))/pi/(1-nu))+B2*(0.125d0*((2-2*nu)*(-2/y1/(1+y2**2/y1
     $     **2)+1/(y1*cos(b)-y3*sin(b))/(1+y2**2/(y1*cos(b)-y3*sin(b))
     $     **2)+((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3
     $     *sin(b))+y2**2*cos(b))+y2**2/(y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))-2*y2**2*(y1
     $     **2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2
     $     **2*cos(b))**2*cos(b))/(1+y2**2*(y1**2+y2**2+y3**2)*sin(b)**2
     $     /(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2)+1/(y1*cos(b)+(y3
     $     +2*a)*sin(b))/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)+((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2
     $     *a)*sin(b))+y2**2*cos(b))+y2**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))-2*y2**2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     /(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*cos(b))/(1
     $     +y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2))+y1*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a))+y1*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y3)*y2-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y3)**2*y2-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*y2-1
     $     /(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*y2)-(y1*cos(b)-y3*sin(b))/(y1**2+y2**2
     $     +y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3
     $     *cos(b))-(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))-y2*(-(y1*cos(b)-y3*sin(b))/(y1**2+y2**2+y3
     $     **2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3
     $     *cos(b))*y2-(y1*cos(b)-y3*sin(b))/(y1**2+y2**2+y3**2)/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*y2-(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y2
     $     -(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *y2))/pi/(1-nu)+1/4*((2-2*nu)*(1-2*nu)*(-1/y1/(1+y2**2/y1**2)
     $     +1/(y1*cos(b)+(y3+2*a)*sin(b))/(1+y2**2/(y1*cos(b)+(y3+2*a)
     $     *sin(b))**2)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1
     $     *(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))+y2**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))
     $     +y2**2*cos(b))-2*y2**2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2
     $     *cos(b))/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1
     $     *cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)**2+(1-2
     $     *nu)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*((-1+2*nu+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+y1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)))-(1-2*nu)*y2**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*((-1+2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*cot(b)+y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+(1-2*nu)*y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*(-a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2
     $     *cot(b)-y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(nu
     $     +a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2-y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1)-(1-2*nu)*cot(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))+(1-2
     $     *nu)*y2**2*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+(1-2*nu)
     $     *y2**2*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)
     $     -a*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)+3*a*y2**2
     $     *(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)+(y3+a)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*((1-2*nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)-a*y1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)))-y2**2*(y3+a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)*((1-2*nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)-a*y1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)))-y2**2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((1-2*nu)
     $     *cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-a
     $     *y1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)))+y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu*y1/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2+a*y1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a))*y2-a*y1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2-1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2))+(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*((-2+2*nu)*cos(b)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     /cos(b))-y2**2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*((-2+2*nu)*cos(b)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     /cos(b))-y2**2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *((-2+2*nu)*cos(b)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))+a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)/cos(b))+y2*(y3+a)*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)*y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)/cos(b))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /cos(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(1.5d0)/cos(b)*y2-2*a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**2/cos(b)*y2))/pi/(1-nu))+B3*(0.125d0*((1-2*nu)*sin(b)
     $     *(1/(y1**2+y2**2+y3**2)**(0.5d0)*y2/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))+1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b)))-2*y2*sin(b)*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     +1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-y2**2*sin(b)*(-1
     $     /(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))*y2-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*y2-1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*y2-1/(y1**2+y2**2+(y3+2*a)**2)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *y2))/pi/(1-nu)+1/4*((1-2*nu)*(-sin(b)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))+y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1+y2/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *y1-(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *y2-(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*y2)-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a))*y1+y1*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(-2*a/(y1**2+y2**2+(y3+2*a)**2)**2*y2-1/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2)+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(sin(b)*(cos(b)-a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))+(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(1+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2))-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)
     $     *sin(b)-a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-(y3+a)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(sin(b)
     $     *a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2-(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(1+a*(y3+2*a)/(y1
     $     **2+y2**2+(y3+2*a)**2))*y2-2*(y1*cos(b)+(y3+2*a)*sin(b))/(y1
     $     **2+y2**2+(y3+2*a)**2)**(5/2)*a*(y3+2*a)*y2+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)*sin(b)-a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*y2+1/(y1**2+y2
     $     **2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(y2**2*cos(b)*sin(b)-a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*y2-1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(2*y2*cos(b)*sin(b)+a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*y2-a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y2)))/pi
     $     /(1-nu))

        e33 = B1*(0.125d0*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)*y3+0.5d0
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)-cos(b)*((1/(y1
     $     **2+y2**2+y3**2)**(0.5d0)*cos(b)*y3-1)/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     -((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y3-((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3
     $     **2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3
     $     *cos(b))**2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y3-cos(b))-(0.5d0
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)+1)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+0.5d0*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(2*y3+4*a)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +cos(b))))/pi/(1-nu)+1/4*((2-2*nu)*((1-2*nu)*(-y2/(y1*cos(b)
     $     +(y3+2*a)*sin(b))**2*sin(b)/(1+y2**2/(y1*cos(b)+(y3+2*a)
     $     *sin(b))**2)+(0.5d0*y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))*(2*y3+4
     $     *a)-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)**2/(y1*(y1
     $     *cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*y1)/(1+y2**2*(y1**2
     $     +y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))
     $     +y2**2*cos(b))**2))*cot(b)-y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+1)-0.5d0*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)+y2*cos(b)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))+0.5d0*y2
     $     *cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a))
     $     +y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*nu/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)+a/(y1**2+y2**2+(y3+2*a)**2))
     $     -0.5d0*y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*nu
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)+a/(y1**2+y2**2
     $     +(y3+2*a)**2))*(2*y3+4*a)+y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(-2*nu/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)-a
     $     /(y1**2+y2**2+(y3+2*a)**2)**2*(2*y3+4*a))+y2*cos(b)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1-2*nu-((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))-a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))
     $     -0.5d0*y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1-2*nu-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))*(2*y3+4*a)-y2*(y3+a)
     $     *cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1-2*nu-((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-a*(y3+2*a)/(y1**2+y2**2
     $     +(y3+2*a)**2))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2
     $     *y3+4*a)+cos(b))+y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(-(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)*(2*y3+4*a)+1)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))
     $     +0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)-a
     $     /(y1**2+y2**2+(y3+2*a)**2)+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**2*(2*y3+4*a)))/pi/(1-nu))+B2*(0.125d0*((-1+2*nu)*sin(b)
     $     *((1/(y1**2+y2**2+y3**2)**(0.5d0)*y3-cos(b))/((y1**2+y2**2+y3
     $     **2)**(0.5d0)-y1*sin(b)-y3*cos(b))-(0.5d0/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-y1*(-1/(y1**2+y2**2
     $     +y3**2)**(1.5d0)*y3+0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *(2*y3+4*a))-sin(b)*((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)
     $     /(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))+(y1*cos(b)-y3*sin(b))*(1/(y1**2+y2**2
     $     +y3**2)**(0.5d0)*cos(b)*y3-1)/(y1**2+y2**2+y3**2)**(0.5d0)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-(y1
     $     *cos(b)-y3*sin(b))*((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)
     $     /(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))*y3-(y1*cos(b)-y3*sin(b))*((y1**2+y2**2
     $     +y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2
     $     +y2**2+y3**2)**(0.5d0)*y3-cos(b))-sin(b)*((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))-(y1*cos(b)+(y3+2*a)*sin(b))*(0.5d0/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)+1)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))+0.5d0*(y1*cos(b)+(y3+2*a)*sin(b))
     $     *((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(2*y3+4*a)+(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b)))/pi/(1-nu)+1/4
     $     *((-2+2*nu)*(1-2*nu)*cot(b)*((0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+1)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)-cos(b)*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2
     $     *y3+4*a)+cos(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b)))+(2-2*nu)*y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+1)+0.5d0*(2-2*nu)*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)+(2-2
     $     *nu)*sin(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     -(2-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(2*y3+4*a)+cos(b))-0.5d0*(2-2*nu)*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *(2*y3+4*a)+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)
     $     *cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-a
     $     *y1/(y1**2+y2**2+(y3+2*a)**2))-0.5d0*(y3+a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)*((1-2*nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)-a*y1/(y1**2+y2**2+(y3+2*a)**2))*(2
     $     *y3+4*a)+(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*nu*y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(0.5d0/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+a*y1/(y1**2+y2**2
     $     +(y3+2*a)**2)**2*(2*y3+4*a))-1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)*sin(b)+((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*cot(b)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(sin(b)-(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))
     $     /(y1**2+y2**2+(y3+2*a)**2)-(y1*cos(b)+(y3+2*a)*sin(b))*((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))))+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)*sin(b)+((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*cot(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*(sin(b)-(y3+2*a)*(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)-(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))))*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))-(y3+a)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *((0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)
     $     +1)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((2-2*nu)
     $     *cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b)))-0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3
     $     +2*a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((2-2*nu)
     $     *cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b)))*(2*y3+4*a)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(
     $     -(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)
     $     +1)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+cos(b)))-0.5d0*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *(sin(b)-(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2
     $     +(y3+2*a)**2)-(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b)))*(2*y3+4*a)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(-(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)-(y3
     $     +2*a)*sin(b)/(y1**2+y2**2+(y3+2*a)**2)+(y3+2*a)*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**2*(2*y3+4*a)
     $     -sin(b)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-(y1*cos(b)+(y3+2*a)
     $     *sin(b))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*(2
     $     *y3+4*a)+1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+0.5d0*(y1
     $     *cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(2*y3+4
     $     *a)+(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+cos(b)))))/pi/(1-nu))+B3*(0.125d0*((2-2*nu)*(y2/(y1
     $     *cos(b)-y3*sin(b))**2*sin(b)/(1+y2**2/(y1*cos(b)-y3*sin(b))
     $     **2)+(y2/(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)
     $     -y3*sin(b))+y2**2*cos(b))*y3+y2*(y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)**2/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2*y1)/(1
     $     +y2**2*(y1**2+y2**2+y3**2)*sin(b)**2/(y1*(y1*cos(b)-y3
     $     *sin(b))+y2**2*cos(b))**2)+y2/(y1*cos(b)+(y3+2*a)*sin(b))**2
     $     *sin(b)/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)-(0.5d0*y2
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))*(2*y3+4*a)-y2*(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2
     $     **2*cos(b))**2*y1)/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     **2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))+y2
     $     *sin(b)*((1/(y1**2+y2**2+y3**2)**(0.5d0)*cos(b)*y3-1)/(y1**2
     $     +y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))-((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)
     $     /(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))*y3-((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)
     $     -y3)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))-(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)*(2*y3+4*a)+1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))+0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(2*y3+4*a)+((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+cos(b))))/pi/(1-nu)+1/4*((2-2*nu)*(-y2
     $     /(y1*cos(b)+(y3+2*a)*sin(b))**2*sin(b)/(1+y2**2/(y1*cos(b)
     $     +(y3+2*a)*sin(b))**2)+(0.5d0*y2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))*(2*y3+4*a)-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2
     $     *y1)/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1
     $     *cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))-(2-2*nu)*y2
     $     *sin(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))
     $     -0.5d0*(2-2*nu)*y2*sin(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(2*y3+4*a)+y2*sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(1+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3
     $     +2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3
     $     +2*a)/(y1**2+y2**2+(y3+2*a)**2))-0.5d0*y2*(y3+a)*sin(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))
     $     *(2*y3+4*a)-y2*(y3+a)*sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))*(0.5d0/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))+y2*(y3+a)*sin(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*((0.5d0/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)+1)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+cos(b))-0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3
     $     +4*a)+a/(y1**2+y2**2+(y3+2*a)**2)-a*(y3+2*a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**2*(2*y3+4*a)))/pi/(1-nu))

        e12 = 0.5d0*B1*(0.125d0*((2-2*nu)*(-2/y1/(1+y2**2/y1**2)+1/(y1
     $     *cos(b)-y3*sin(b))/(1+y2**2/(y1*cos(b)-y3*sin(b))**2)+((y1**2
     $     +y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2
     $     *cos(b))+y2**2/(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1
     $     *cos(b)-y3*sin(b))+y2**2*cos(b))-2*y2**2*(y1**2+y2**2+y3**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2
     $     *cos(b))/(1+y2**2*(y1**2+y2**2+y3**2)*sin(b)**2/(y1*(y1
     $     *cos(b)-y3*sin(b))+y2**2*cos(b))**2)+1/(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)+((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))+y2**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))-2*y2**2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     /(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*cos(b))/(1
     $     +y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2))-y1*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a))-y1*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y3)*y2-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y3)**2*y2-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*y2-1
     $     /(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*y2)-cos(b)*(((y1**2+y2**2+y3**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2
     $     **2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))+((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b)))-y2*cos(b)*(1/(y1**2+y2**2+y3**2)*sin(b)*y2/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-((y1**2+y2**2
     $     +y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)**(1.5d0)/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))*y2-((y1**2+y2
     $     **2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)/((y1**2+y2
     $     **2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*y2+1/(y1**2+y2**2
     $     +(y3+2*a)**2)*sin(b)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y2
     $     -((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*y2))/pi/(1-nu)+1/4*((-2+2*nu)*(1-2*nu)*(
     $     -1/y1/(1+y2**2/y1**2)+1/(y1*cos(b)+(y3+2*a)*sin(b))/(1+y2**2
     $     /(y1*cos(b)+(y3+2*a)*sin(b))**2)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))+y2**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1
     $     *(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))-2*y2**2*(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))**2*cos(b))/(1+y2**2*(y1**2+y2**2+(y3+2
     $     *a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))**2))*cot(b)**2+(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*((1-2*nu-a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*cot(b)-y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)))-(1-2*nu)*y2**2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((1-2*nu-a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)-y1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+(1-2*nu)
     $     *y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)*y2*cot(b)+y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)**2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2+y2/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*y1)+(1-2*nu)*cos(b)*cot(b)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))-(1-2*nu)*y2**2*cos(b)*cot(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-(1-2*nu)*y2**2*cos(b)*cot(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)+a*(y3+a)*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-3*a*y2**2*(y3+a)*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(5/2)+(y3+a)/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *((-1+2*nu)*cot(b)+y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*y1/(y1**2
     $     +y2**2+(y3+2*a)**2))-y2**2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*((-1+2
     $     *nu)*cot(b)+y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2
     $     *nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*y1/(y1**2+y2**2
     $     +(y3+2*a)**2))-y2**2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2*nu)*cot(b)+y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))+a*y1/(y1**2+y2**2+(y3+2*a)**2))
     $     +y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-y2/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*y1-2*a*y1/(y1**2+y2**2+(y3+2*a)**2)**2*y2)+(y3
     $     +a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)
     $     -a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))-a*(y3+2*a)
     $     *cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2))-y2**2*(y3+a)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))-a*(y3+2*a)
     $     *cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2))-y2**2*(y3+a)/(y1**2
     $     +y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(cos(b)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))-a*(y3+2*a)*cos(b)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2))+y2*(y3+a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(-cos(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2+cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)*y2*((1-2*nu)*cos(b)-a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))*cot(b)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2
     $     *cot(b)+(2-2*nu)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)*y2
     $     *cos(b))+2*a*(y3+2*a)*cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **2*y2))/pi/(1-nu))+0.5d0*B2*(0.125d0*((-1+2*nu)*(1/(y1**2+y2
     $     **2+y3**2)**(0.5d0)*y2/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)-cos(b)*(1/(y1**2+y2**2+y3**2)**(0.5d0)
     $     *y2/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))+1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))))+y1**2*(-1/(y1**2+y2**2
     $     +y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)*y2-1/(y1
     $     **2+y2**2+y3**2)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)**2*y2-1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*y2-1/(y1**2+y2**2+(y3+2*a)**2)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*y2)+(y1*cos(b)-y3
     $     *sin(b))/(y1**2+y2**2+y3**2)*sin(b)*y2/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))-(y1*cos(b)-y3*sin(b))*((y1**2
     $     +y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y2-(y1*cos(b)-y3*sin(b))*((y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)-y1)/(y1**2+y2**2+y3**2)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))**2*y2+(y1*cos(b)+(y3+2*a)*sin(b))/(y1
     $     **2+y2**2+(y3+2*a)**2)*sin(b)*y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y2-(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*y2)/pi/(1-nu)+1/4*((1-2*nu)*(((2
     $     -2*nu)*cot(b)**2+nu)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-((2-2*nu)*cot(b)
     $     **2+1)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-(1-2
     $     *nu)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2
     $     *nu)*y1*cot(b)+nu*(y3+2*a)-a+a*y1*cot(b)/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y2+(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*(-a*y1*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y2-y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)**2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y2-y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2)+(1
     $     -2*nu)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*((y1*cos(b)+(y3+2*a)*sin(b))*cos(b)-a
     $     *((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/cos(b))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2-(1-2*nu)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(-a/(y1**2+y2**2+(y3+2
     $     *a)**2)*sin(b)*y2/cos(b)+a*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /cos(b)*y2)+3*a*y2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(5/2)*y1-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*(2*nu+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*y1
     $     *cot(b)+a)-y1**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))-a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2+(y3+a)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)*((1-2*nu)*y1*cot(b)+a)*y2+y1**2/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*y2+y1**2
     $     /(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*y2+y1**2/(y1**2+y2**2+(y3+2*a)**2)**2/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a*y2+3*a*y1**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(5/2)*y2)-(y3+a)*cot(b)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(-cos(b)*sin(b)+a
     $     *y1*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)+((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2+(y3
     $     +a)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))*(-3*a*y1*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(5
     $     /2)/cos(b)*y2+1/(y1**2+y2**2+(y3+2*a)**2)*sin(b)*y2*((2-2*nu)
     $     *cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b)))
     $     -((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)/cos(b)))*y2+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(-1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*y2/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))+((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/cos(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *y2+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)*y2)))/pi/(1-nu))
     $     +0.5d0*B3*(0.125d0*sin(b)*(((y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)-y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))
     $     /pi/(1-nu)+0.125d0*y2*sin(b)*(1/(y1**2+y2**2+y3**2)*sin(b)*y2
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-((y1**2
     $     +y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y2-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3
     $     **2)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*y2
     $     +1/(y1**2+y2**2+(y3+2*a)**2)*sin(b)*y2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*y2-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)
     $     /(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*y2)/pi/(1-nu)+1/4*((1
     $     -2*nu)*(1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(1+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-y2**2/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y2**2/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)-cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+y2**2*cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y2**2
     $     *cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0))-(y3+a)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)
     $     +1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a))+y2**2*(y3+a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(a/(y1**2+y2**2+(y3+2*a)
     $     **2)+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a))-y2*(y3+a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(-2*a/(y1**2+y2**2+(y3+2
     $     *a)**2)**2*y2-1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2)+(y3+a)*cos(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))
     $     -y2**2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1
     $     **2+y2**2+(y3+2*a)**2))-y2**2*(y3+a)*cos(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3
     $     +2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3
     $     +2*a)/(y1**2+y2**2+(y3+2*a)**2))+y2*(y3+a)*cos(b)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y2-2*a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**2*y2))
     $     /pi/(1-nu))+0.5d0*B1*(0.125d0*((1-2*nu)*(1/(y1**2+y2**2+y3
     $     **2)**(0.5d0)*y1/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)-cos(b)*((1/(y1**2+y2**2+y3**2)**(0.5d0)*y1
     $     -sin(b))/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     +(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))))-y2**2
     $     *(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y3)*y1-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y3)**2*y1-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*y1-1/(y1**2+y2**2+(y3
     $     +2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*y1
     $     -cos(b)*(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))*y1-1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     **2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))-1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*y1-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b)))))/pi/(1-nu)+1/4*((1-2*nu)*(((2-2*nu)*cot(b)**2-nu)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)-((2-2*nu)*cot(b)**2+1-2*nu)*cos(b)*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))+(1-2*nu)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(y1*cot(b)*(1
     $     -2*nu-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+nu*(y3+2*a)-a+y2
     $     **2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1-(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)*((1-2*nu-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+a
     $     *y1**2*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-y2**2/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(nu+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -y2**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0)*y1)-(1-2*nu)*cos(b)*cot(b)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+(1-2*nu)*(y1
     $     *cos(b)+(y3+2*a)*sin(b))*cot(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1-sin(b))+(1-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))
     $     *cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1-a*(y3+a)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)+3*a*y1**2*(y3+a)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)-(y3+a)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(-2*nu+1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*((1-2*nu)*y1*cot(b)-a)+y2**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*y2**2/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((1-2*nu)*y1
     $     *cot(b)-a)*y1+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1-2*nu)
     $     *cot(b)-y2**2/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))*y1-y2**2/(y1**2+y2**2+(y3+2*a)**2)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))*y1-y2**2/(y1**2+y2**2+(y3+2*a)**2)**2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a*y1-3*a*y2**2
     $     /(y1**2+y2**2+(y3+2*a)**2)**(5/2)*y1)-(y3+a)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)
     $     **2-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*(y1*cos(b)
     $     +(y3+2*a)*sin(b))*cot(b)+a*cos(b))+a*(y3+2*a)*(y1*cos(b)+(y3
     $     +2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)**2-a*(y1
     $     *cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+(y3+a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((1-2*nu)*(y1
     $     *cos(b)+(y3+2*a)*sin(b))*cot(b)+a*cos(b))*y1-1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(1-2*nu)*cos(b)*cot(b)+a*(y3+2*a)
     $     *cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-3*a*(y3+2
     $     *a)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(5/2)*y1+1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2
     $     **2*cos(b)**2-a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a))*y1+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(y2**2*cos(b)**2-a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a))*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1-sin(b))-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(-a*cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)+a*(y1
     $     *cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     *y1-a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2
     $     *a)**2)*cos(b)*y1)))/pi/(1-nu))+0.5d0*B2*(0.125d0*((2-2*nu)
     $     *(2*y2/y1**2/(1+y2**2/y1**2)-y2/(y1*cos(b)-y3*sin(b))**2
     $     *cos(b)/(1+y2**2/(y1*cos(b)-y3*sin(b))**2)+(y2/(y1**2+y2**2
     $     +y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2
     $     *cos(b))*y1-y2*(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1
     $     *cos(b)-y3*sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)-y3*sin(b)))
     $     /(1+y2**2*(y1**2+y2**2+y3**2)*sin(b)**2/(y1*(y1*cos(b)-y3
     $     *sin(b))+y2**2*cos(b))**2)-y2/(y1*cos(b)+(y3+2*a)*sin(b))**2
     $     *cos(b)/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)+(y2/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))*y1-y2*(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))**2*(2*y1*cos(b)+(y3+2*a)*sin(b)))/(1+y2**2*(y1**2+y2
     $     **2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2
     $     **2*cos(b))**2))+y2*(1/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y3)+1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a))+y1*y2
     $     *(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y3)*y1-1/(y1**2+y2**2+y3**2)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y3)**2*y1-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*y1-1/(y1**2+y2**2+(y3
     $     +2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*y1)
     $     -y2*(cos(b)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))-(y1*cos(b)-y3*sin(b))/(y1**2
     $     +y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))*y1-(y1*cos(b)-y3*sin(b))/(y1**2+y2**2+y3
     $     **2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3
     $     *cos(b))**2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))+cos(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1-(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))))/pi/(1-nu)+1
     $     /4*((2-2*nu)*(1-2*nu)*(y2/y1**2/(1+y2**2/y1**2)-y2/(y1*cos(b)
     $     +(y3+2*a)*sin(b))**2*cos(b)/(1+y2**2/(y1*cos(b)+(y3+2*a)
     $     *sin(b))**2)+(y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     /(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))*y1-y2*(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)+(y3+2*a)*sin(b)))/(1
     $     +y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)**2-(1-2*nu)*y2/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2*nu+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+y1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1+(1-2*nu)*y2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*y1*cot(b)+1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     -y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(nu+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0))+(1-2*nu)*y2*cot(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))
     $     *(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+(1-2*nu)*y2
     $     *cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)*y1+3
     $     *a*y2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)*y1-y2
     $     *(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)*((1-2*nu)*cot(b)-2*nu*y1/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-a*y1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)))*y1-y2*(y3+a)/(y1**2
     $     +y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)**2*((1-2*nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)-a*y1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)))*y1+y2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-2
     $     *nu/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)+2*nu*y1**2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a))+a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a))-a*y1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(-1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1-1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1))-y2*(y3+a)*cot(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*((-2+2*nu)*cos(b)+((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)/cos(b))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     /cos(b))*y1-y2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*((-2+2*nu)*cos(b)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     /cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+y2
     $     *(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*y1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/cos(b))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)*y1-2*a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**2/cos(b)*y1))/pi/(1-nu))
     $     +0.5d0*B3*(0.125d0*((1-2*nu)*sin(b)*((1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y1-sin(b))/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)
     $     -y3*cos(b))+(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b)))-y2**2*sin(b)*(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))*y1-1/(y1**2+y2
     $     **2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)
     $     -y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))-1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1-1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y1-sin(b))))/pi/(1-nu)+1/4*((1-2*nu)*(-sin(b)*(1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-1/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(1+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))+y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)+cos(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-(y1
     $     *cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))-(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(1.5d0)*y1)+(y3+a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a))-y1**2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a))+y1*(y3+a)/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(-2*a/(y1**2+y2**2+(y3+2*a)**2)**2*y1-1/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*y1)+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(sin(b)*(cos(b)-a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0))+(y1*cos(b)+(y3+2*a)*sin(b))
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1+a*(y3+2*a)/(y1**2+y2
     $     **2+(y3+2*a)**2))-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(y2**2*cos(b)*sin(b)-a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(sin(b)*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y1+cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1+a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))-(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(1+a*(y3+2*a)/(y1
     $     **2+y2**2+(y3+2*a)**2))*y1-2*(y1*cos(b)+(y3+2*a)*sin(b))/(y1
     $     **2+y2**2+(y3+2*a)**2)**(5/2)*a*(y3+2*a)*y1+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)*sin(b)-a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*y1+1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(y2**2*cos(b)*sin(b)-a*(y1
     $     *cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*(1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))-1/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(-a*cos(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     +a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     *y1-a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)
     $     *cos(b)*y1)))/pi/(1-nu))

        e13 = 0.5d0*B1*(0.125d0*((2-2*nu)*(y2/(y1*cos(b)-y3*sin(b))**2
     $     *sin(b)/(1+y2**2/(y1*cos(b)-y3*sin(b))**2)+(y2/(y1**2+y2**2
     $     +y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2
     $     *cos(b))*y3+y2*(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)**2/(y1*(y1
     $     *cos(b)-y3*sin(b))+y2**2*cos(b))**2*y1)/(1+y2**2*(y1**2+y2**2
     $     +y3**2)*sin(b)**2/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2)
     $     -y2/(y1*cos(b)+(y3+2*a)*sin(b))**2*sin(b)/(1+y2**2/(y1*cos(b)
     $     +(y3+2*a)*sin(b))**2)+(0.5d0*y2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))*(2*y3+4*a)-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2
     $     *y1)/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1
     $     *cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))-y1*y2*(-1/(y1**2
     $     +y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)*y3-1
     $     /(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y3)**2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y3-1)-0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)*(2*y3+4*a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1))-y2*cos(b)*(1/(y1**2
     $     +y2**2+y3**2)*sin(b)*y3/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)-y1)
     $     /(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))*y3-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)
     $     -y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))+0.5d0/(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     *(2*y3+4*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))-0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(2*y3+4*a)
     $     -((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+cos(b))))/pi/(1-nu)+1/4*((-2+2*nu)*(1-2
     $     *nu)*(-y2/(y1*cos(b)+(y3+2*a)*sin(b))**2*sin(b)/(1+y2**2/(y1
     $     *cos(b)+(y3+2*a)*sin(b))**2)+(0.5d0*y2/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))*(2*y3+4*a)-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2
     $     *y1)/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1
     $     *cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)**2-(1-2
     $     *nu)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((1-2
     $     *nu-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)-y1/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+1)+(1-2*nu)*y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*(0.5d0*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *(2*y3+4*a)*cot(b)+y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)**2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+0.5d0*y1/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(1.5d0)*(2*y3+4*a))-(1-2*nu)*y2*cos(b)*cot(b)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))-0.5d0*(1-2
     $     *nu)*y2*cos(b)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *(2*y3+4*a)+a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2*cot(b)
     $     -1.5d0*a*y2*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)*(2
     $     *y3+4*a)+y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*((-1+2*nu)*cot(b)+y1/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))+a*y1/(y1**2+y2**2+(y3+2*a)**2))-0.5d0*y2
     $     *(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)*((-1+2*nu)*cot(b)+y1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))+a*y1/(y1**2+y2**2+(y3+2*a)**2))*(2*y3+4*a)-y2
     $     *(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2*nu)*cot(b)+y1/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))+a*y1/(y1**2+y2**2+(y3+2*a)**2))*(0.5d0/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+y2*(y3+a)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*(-y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)-0.5d0*y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)-a*y1/(y1**2+y2**2+(y3+2*a)
     $     **2)**2*(2*y3+4*a))+y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*cot(b)+(2-2*nu)*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)*cos(b))-a*(y3+2*a)*cos(b)*cot(b)/(y1**2
     $     +y2**2+(y3+2*a)**2))-0.5d0*y2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*sin(b)-y1)*cos(b))-a*(y3+2*a)*cos(b)*cot(b)/(y1
     $     **2+y2**2+(y3+2*a)**2))*(2*y3+4*a)-y2*(y3+a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(cos(b)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))-a*(y3+2*a)*cos(b)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2))*(0.5d0/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))+y2*(y3+a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(-cos(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)*((1-2*nu)*cos(b)-a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))*cot(b)+(2-2*nu)*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*sin(b)-y1)*cos(b))*(0.5d0/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))+cos(b)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*((0.5d0/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)+1)*((1-2
     $     *nu)*cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)
     $     +0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)*cot(b)+0.5d0
     $     *(2-2*nu)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)*(2*y3+4
     $     *a)*cos(b))-a*cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)+a*(y3+2
     $     *a)*cos(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**2*(2*y3+4*a)))
     $     /pi/(1-nu))+0.5d0*B2*(0.125d0*((-1+2*nu)*((1/(y1**2+y2**2+y3
     $     **2)**(0.5d0)*y3-1)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+(0.5d0
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-cos(b)*((1/(y1**2+y2**2+y3
     $     **2)**(0.5d0)*y3-cos(b))/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))+(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+cos(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))))+y1**2*(-1/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)*y3-1/(y1**2+y2**2
     $     +y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)**2*(1/(y1
     $     **2+y2**2+y3**2)**(0.5d0)*y3-1)-0.5d0/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2
     $     *y3+4*a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(0.5d0/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(2*y3+4*a)+1))-sin(b)*((y1**2+y2**2+y3**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2
     $     **2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))+(y1*cos(b)-y3
     $     *sin(b))/(y1**2+y2**2+y3**2)*sin(b)*y3/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))-(y1*cos(b)-y3*sin(b))*((y1**2
     $     +y2**2+y3**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y3-(y1*cos(b)-y3*sin(b))*((y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)-y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))+sin(b)*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     +0.5d0*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)
     $     *sin(b)*(2*y3+4*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))-0.5d0*(y1*cos(b)+(y3+2*a)*sin(b))
     $     *((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(2*y3+4*a)-(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b)))/pi/(1-nu)+1/4*((1
     $     -2*nu)*(((2-2*nu)*cot(b)**2+nu)*(0.5d0/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(2*y3+4*a)+1)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)-((2-2*nu)*cot(b)**2+1)*cos(b)*(0.5d0/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-(1-2*nu)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((-1+2*nu)*y1
     $     *cot(b)+nu*(y3+2*a)-a+a*y1*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y1**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)))*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+(1-2*nu)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu-0.5d0*a*y1*cot(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)-y1**2/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(nu+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+1)-0.5d0*y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3
     $     +4*a))+(1-2*nu)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*((y1*cos(b)+(y3+2*a)*sin(b))
     $     *cos(b)-a*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))*(0.5d0/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))-(1-2*nu)*cot(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)*sin(b)-0.5d0*a/(y1**2+y2**2+(y3+2*a)**2)
     $     *sin(b)*(2*y3+4*a)/cos(b)+0.5d0*a*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /cos(b)*(2*y3+4*a))-a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1
     $     *cot(b)+1.5d0*a*y1*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(5/2)*(2*y3+4*a)+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*(2*nu+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*y1
     $     *cot(b)+a)-y1**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))-a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0))-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*(2*nu+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*y1
     $     *cot(b)+a)-y1**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))-a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+1)+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(
     $     -0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((1-2*nu)*y1*cot(b)
     $     +a)*(2*y3+4*a)+0.5d0*y1**2/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))*(2*y3+4*a)+y1**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+0.5d0*y1**2
     $     /(y1**2+y2**2+(y3+2*a)**2)**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a*(2*y3+4*a)+1.5d0*a*y1**2/(y1**2+y2**2+(y3
     $     +2*a)**2)**(5/2)*(2*y3+4*a))+cot(b)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(-cos(b)*sin(b)+a*y1
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)+((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))))-(y3+a)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(-cos(b)*sin(b)+a*y1
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)+((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+cos(b))+(y3+a)*cot(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)/cos(b)*y1-1.5d0*a*y1*(y3+2*a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(5/2)/cos(b)*(2*y3+4*a)+0.5d0/(y1**2+y2**2+(y3+2
     $     *a)**2)*sin(b)*(2*y3+4*a)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)/cos(b)))-0.5d0*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((2-2
     $     *nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b)))*(2
     $     *y3+4*a)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(-(0.5d0/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)+1)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/cos(b))+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+cos(b))+0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)
     $     *(2*y3+4*a))))/pi/(1-nu))+0.5d0*B3*(0.125d0*y2*sin(b)*(1/(y1
     $     **2+y2**2+y3**2)*sin(b)*y3/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)-y1)
     $     /(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))*y3-((y1**2+y2**2+y3**2)**(0.5d0)*sin(b)
     $     -y1)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))+0.5d0/(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     *(2*y3+4*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))-0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)-y1)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(2*y3+4*a)
     $     -((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)-y1)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+cos(b)))/pi/(1-nu)+1/4*((1-2*nu)*(-y2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(1+a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(2*y3+4*a)+1)-0.5d0*y2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2
     $     *y3+4*a)+y2*cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2
     $     *y3+4*a)+cos(b))+0.5d0*y2*cos(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*(2*y3+4*a))-y2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a))+0.5d0*y2*(y3+a)/(y1**2+y2**2+(y3+2
     $     *a)**2)**(1.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a))*(2*y3+4*a)-y2*(y3+a)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(-a/(y1**2+y2**2+(y3+2*a)**2)**2
     $     *(2*y3+4*a)-1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1))+y2
     $     *cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2))-0.5d0*y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3
     $     +2*a)/(y1**2+y2**2+(y3+2*a)**2))*(2*y3+4*a)-y2*(y3+a)*cos(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +cos(b))+y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*((0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*(2
     $     *y3+4*a)+1)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     -((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))-0.5d0*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     -a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**2*(2*y3+4*a)))/pi/(1
     $     -nu))+0.5d0*B1*(0.125d0*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)
     $     *y1+1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1-cos(b)*(1/(y1**2
     $     +y2**2+y3**2)*cos(b)*y1/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))-((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)
     $     /(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))*y1-((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)
     $     -y3)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y1-sin(b))-1/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1+((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))))/pi/(1-nu)+1/4*((2-2*nu)*((1-2*nu)*(y2/y1**2/(1+y2
     $     **2/y1**2)-y2/(y1*cos(b)+(y3+2*a)*sin(b))**2*cos(b)/(1+y2**2
     $     /(y1*cos(b)+(y3+2*a)*sin(b))**2)+(y2/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))*y1-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1
     $     *(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)
     $     +(y3+2*a)*sin(b)))/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     **2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)
     $     -y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(2*nu+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2-y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1+y2*cos(b)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+y2*cos(b)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)*y1)-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*(2*nu/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)+a/(y1**2+y2**2+(y3+2*a)**2))*y1+y2*(y3+a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(-2*nu/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-2
     $     *a/(y1**2+y2**2+(y3+2*a)**2)**2*y1)-y2*(y3+a)*cos(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1-2*nu-((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))-a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2))*y1-y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(1-2*nu-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))-a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))*(1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))+y2*(y3+a)*cos(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(-1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))
     $     +((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)*y1+2*a*(y3+2*a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**2*y1))/pi/(1-nu))+0.5d0*B2*(0.125d0*((-1+2*nu)
     $     *sin(b)*((1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*(-1/(y1**2
     $     +y2**2+y3**2)**(1.5d0)*y1+1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y1)+cos(b)*((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)
     $     /(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))+(y1*cos(b)-y3*sin(b))/(y1**2+y2**2+y3
     $     **2)*cos(b)*y1/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3
     $     *cos(b))-(y1*cos(b)-y3*sin(b))*((y1**2+y2**2+y3**2)**(0.5d0)
     $     *cos(b)-y3)/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))*y1-(y1*cos(b)-y3*sin(b))*((y1
     $     **2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     **2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y1-sin(b))-cos(b)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))-(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2
     $     +y2**2+(y3+2*a)**2)*cos(b)*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1+(y1*cos(b)+(y3+2
     $     *a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*y1-sin(b)))/pi/(1-nu)+1/4*((-2+2*nu)
     $     *(1-2*nu)*cot(b)*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-cos(b)*(1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))-(2-2*nu)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0))+(2-2*nu)*y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+(2-2*nu)*y1**2
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)+(2-2*nu)*cos(b)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))-(2-2*nu)*(y1*cos(b)+(y3+2*a)
     $     *sin(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     *(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))-(2-2*nu)
     $     *(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*y1-(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *((1-2*nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)-a*y1/(y1**2+y2**2+(y3+2*a)**2))*y1+(y3+a)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(-2*nu/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)+2*nu*y1**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-a/(y1
     $     **2+y2**2+(y3+2*a)**2)+2*a*y1**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **2)+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))**2*(cos(b)*sin(b)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b)))+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(sin(b)-(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2
     $     +y2**2+(y3+2*a)**2)-(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y1
     $     *cot(b)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b)))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((2-2*nu)
     $     *cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b)))*y1+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(-1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*cos(b)*y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*y1-sin(b)))-a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(sin(b)-(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2
     $     +y2**2+(y3+2*a)**2)-(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b)))*y1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(
     $     -(y3+2*a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)+2*(y3+2*a)*(y1
     $     *cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**2*y1
     $     -cos(b)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))-(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+(y1*cos(b)
     $     +(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1+(y1*cos(b)
     $     +(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*y1-sin(b)))))/pi/(1-nu))+0.5d0*B3
     $     *(0.125d0*((2-2*nu)*(-y2/(y1*cos(b)-y3*sin(b))**2*cos(b)/(1
     $     +y2**2/(y1*cos(b)-y3*sin(b))**2)+(y2/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))*y1
     $     -y2*(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3
     $     *sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)-y3*sin(b)))/(1+y2**2
     $     *(y1**2+y2**2+y3**2)*sin(b)**2/(y1*(y1*cos(b)-y3*sin(b))+y2
     $     **2*cos(b))**2)+y2/(y1*cos(b)+(y3+2*a)*sin(b))**2*cos(b)/(1
     $     +y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)-(y2/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))*y1-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1
     $     *(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)
     $     +(y3+2*a)*sin(b)))/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     **2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))+y2
     $     *sin(b)*(1/(y1**2+y2**2+y3**2)*cos(b)*y1/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))-((y1**2+y2**2+y3**2)**(0.5d0)
     $     *cos(b)-y3)/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))*y1-((y1**2+y2**2+y3**2)
     $     **(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2
     $     **2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2
     $     +y3**2)**(0.5d0)*y1-sin(b))-1/(y1**2+y2**2+(y3+2*a)**2)
     $     *cos(b)*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y1+((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y1
     $     -sin(b))))/pi/(1-nu)+1/4*((2-2*nu)*(y2/y1**2/(1+y2**2/y1**2)
     $     -y2/(y1*cos(b)+(y3+2*a)*sin(b))**2*cos(b)/(1+y2**2/(y1*cos(b)
     $     +(y3+2*a)*sin(b))**2)+(y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))*y1-y2
     $     *(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2*(2*y1*cos(b)+(y3+2*a)*sin(b)))
     $     /(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)
     $     +(y3+2*a)*sin(b))+y2**2*cos(b))**2))-(2-2*nu)*y2*sin(b)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*y1-sin(b))-(2-2*nu)*y2*sin(b)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1-y2*(y3+a)*sin(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))
     $     *y1-y2*(y3+a)*sin(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(1+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1
     $     **2+y2**2+(y3+2*a)**2))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *y1-sin(b))+y2*(y3+a)*sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)
     $     +a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*y1-sin(b))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y1-2*a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**2*y1))
     $     /pi/(1-nu))

        e23 = 0.5d0*B1*(0.125d0*((1-2*nu)*((1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-1)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)+(0.5d0/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)-cos(b)*((1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)
     $     -y3*cos(b))+(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+cos(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b))))-y2**2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)/((y1
     $     **2+y2**2+y3**2)**(0.5d0)-y3)*y3-1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)**2*(1/(y1**2+y2
     $     **2+y3**2)**(0.5d0)*y3-1)-0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*y3+4
     $     *a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+1)-cos(b)*(-1/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y3-1/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))-0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(2*y3+4*a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +cos(b)))))/pi/(1-nu)+1/4*((1-2*nu)*(((2-2*nu)*cot(b)**2-nu)
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-((2-2*nu)*cot(b)**2+1
     $     -2*nu)*cos(b)*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3
     $     +4*a)+cos(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b)))+(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(y1*cot(b)*(1-2*nu-a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))+nu*(y3+2*a)-a+y2**2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+1)-(1-2*nu)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *(0.5d0*a*y1*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3
     $     +4*a)+nu-y2**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2
     $     *(nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)-0.5d0*y2**2/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*(2*y3+4*a))-(1-2*nu)*sin(b)*cot(b)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)
     $     +a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+(1-2*nu)*(y1*cos(b)
     $     +(y3+2*a)*sin(b))*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+cos(b))+0.5d0*(1-2*nu)*(y1*cos(b)+(y3+2*a)
     $     *sin(b))*cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3
     $     +4*a)-a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y1*cot(b)+1.5d0*a
     $     *y1*(y3+a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)*(2*y3+4*a)
     $     +1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-2*nu+1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*y1*cot(b)-a)+y2**2
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+a*y2**2/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0))-(y3+a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(-2*nu+1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*y1*cot(b)-a)+y2**2
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+a*y2**2/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0))*(0.5d0
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+(y3+a)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(-0.5d0/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*((1-2*nu)*y1*cot(b)-a)*(2*y3+4*a)
     $     -0.5d0*y2**2/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))*(2*y3+4*a)-y2**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(2
     $     *nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*(0.5d0/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)-0.5d0*y2**2/(y1**2+y2**2
     $     +(y3+2*a)**2)**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     *a*(2*y3+4*a)-1.5d0*a*y2**2/(y1**2+y2**2+(y3+2*a)**2)**(5/2)
     $     *(2*y3+4*a))+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(cos(b)**2-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((1-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)+a
     $     *cos(b))+a*(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0)-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(y2**2*cos(b)**2-a*(y1*cos(b)+(y3+2*a)*sin(b))
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a)))-(y3+a)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)**2-1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((1-2*nu)*(y1*cos(b)+(y3
     $     +2*a)*sin(b))*cot(b)+a*cos(b))+a*(y3+2*a)*(y1*cos(b)+(y3+2*a)
     $     *sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)**2-a*(y1
     $     *cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +cos(b))+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *((1-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)+a*cos(b))*(2*y3
     $     +4*a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1-2*nu)*sin(b)
     $     *cot(b)+a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)+a*(y3+2*a)*sin(b)*cot(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)-1.5d0*a*(y3+2*a)*(y1*cos(b)+(y3+2*a)
     $     *sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(5/2)*(2*y3+4*a)
     $     +0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)**2
     $     -a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a))*(2*y3+4*a)+1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2
     $     *(y2**2*cos(b)**2-a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+cos(b))-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(-a*sin(b)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     +0.5d0*a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(1.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)*(2*y3+4*a)-a*(y1*cos(b)+(y3+2*a)*sin(b))*cot(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(0.5d0/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*cos(b)*(2*y3+4*a)+1))))/pi/(1-nu))+0.5d0*B2
     $     *(0.125d0*((2-2*nu)*(y2/(y1*cos(b)-y3*sin(b))**2*sin(b)/(1+y2
     $     **2/(y1*cos(b)-y3*sin(b))**2)+(y2/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))*y3
     $     +y2*(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)**2/(y1*(y1*cos(b)-y3
     $     *sin(b))+y2**2*cos(b))**2*y1)/(1+y2**2*(y1**2+y2**2+y3**2)
     $     *sin(b)**2/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))**2)-y2/(y1
     $     *cos(b)+(y3+2*a)*sin(b))**2*sin(b)/(1+y2**2/(y1*cos(b)+(y3+2
     $     *a)*sin(b))**2)+(0.5d0*y2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))*(2*y3+4
     $     *a)-y2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)**2/(y1*(y1
     $     *cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*y1)/(1+y2**2*(y1**2
     $     +y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))
     $     +y2**2*cos(b))**2))+y1*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y3)*y3-1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y3)**2*(1/(y1**2+y2
     $     **2+y3**2)**(0.5d0)*y3-1)-0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(2*y3+4
     $     *a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+1))-y2*(-sin(b)/(y1**2+y2**2+y3**2)
     $     **(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     -(y1*cos(b)-y3*sin(b))/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))*y3-(y1*cos(b)-y3
     $     *sin(b))/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))+sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))-0.5d0*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(2*y3+4*a)-(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))))/pi/(1
     $     -nu)+1/4*((2-2*nu)*(1-2*nu)*(-y2/(y1*cos(b)+(y3+2*a)*sin(b))
     $     **2*sin(b)/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)+(0.5d0*y2
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))*(2*y3+4*a)-y2*(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2
     $     **2*cos(b))**2*y1)/(1+y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)
     $     **2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)
     $     **2-(1-2*nu)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*((-1+2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))*cot(b)+y1
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*(nu+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+1)+(1-2*nu)*y2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)*(-0.5d0*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(2*y3+4*a)*cot(b)-y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)**2*(nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+1)-0.5d0*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a))+(1-2*nu)*y2
     $     *cot(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /cos(b))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +cos(b))+0.5d0*(1-2*nu)*y2*cot(b)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)/cos(b)*(2*y3+4*a)-a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y2*cot(b)+1.5d0*a*y2*(y3+a)*cot(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(5/2)*(2*y3+4*a)+y2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*((1-2
     $     *nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)-a*y1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     +y3+2*a)))-0.5d0*y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*((1-2*nu)*cot(b)
     $     -2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-a*y1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)))*(2
     $     *y3+4*a)-y2*(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*((1-2*nu)*cot(b)-2*nu
     $     *y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)-a*y1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)))
     $     *(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+1)+y2
     $     *(y3+a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)+y3+2*a)*(2*nu*y1/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+1)+0.5d0*a*y1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+1/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)+y3+2*a))*(2*y3+4*a)-a*y1/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(-0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(2*y3+4*a)-1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3
     $     +2*a)**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +1)))+y2*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*((-2+2
     $     *nu)*cos(b)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))+a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)/cos(b))-0.5d0*y2*(y3+a)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*((-2+2*nu)
     $     *cos(b)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))+a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)/cos(b))*(2*y3+4*a)-y2*(y3
     $     +a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*((-2+2
     $     *nu)*cos(b)+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/cos(b))+a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)/cos(b))*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))+y2*(y3+a)*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*((0.5d0/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)*(2*y3+4*a)+1)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/cos(b))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)/cos(b))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(2*y3+4*a)+cos(b))-0.5d0*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/cos(b)
     $     *(2*y3+4*a)+a/(y1**2+y2**2+(y3+2*a)**2)/cos(b)-a*(y3+2*a)/(y1
     $     **2+y2**2+(y3+2*a)**2)**2/cos(b)*(2*y3+4*a)))/pi/(1-nu))
     $     +0.5d0*B3*(0.125d0*((1-2*nu)*sin(b)*((1/(y1**2+y2**2+y3**2)
     $     **(0.5d0)*y3-cos(b))/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)
     $     -y3*cos(b))+(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+cos(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3
     $     +2*a)*cos(b)))-y2**2*sin(b)*(-1/(y1**2+y2**2+y3**2)**(1.5d0)
     $     /((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))*y3-1/(y1
     $     **2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1
     $     *sin(b)-y3*cos(b))**2*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y3
     $     -cos(b))-0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(2*y3+4
     $     *a)-1/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(0.5d0/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))))/pi/(1-nu)+1
     $     /4*((1-2*nu)*(-sin(b)*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(2*y3+4*a)+cos(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+y1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)**2*(1+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+1)+0.5d0*y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a)+sin(b)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))
     $     *(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-(y1*cos(b)+(y3
     $     +2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4
     $     *a)+cos(b))-0.5d0*(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(1.5d0)*(2*y3+4*a))+y1/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a))-0.5d0*y1*(y3+a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*(a/(y1**2+y2**2+(y3+2*a)**2)+1/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a))*(2*y3+4*a)+y1*(y3+a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(-a/(y1**2+y2**2+(y3+2*a)
     $     **2)**2*(2*y3+4*a)-1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)**2*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +1))-1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(sin(b)*(cos(b)-a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))+(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(1+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))-1/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)*sin(b)-a
     $     *(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)))+(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2*(sin(b)*(cos(b)-a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))+(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*(1+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))-1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(y2**2*cos(b)*sin(b)
     $     -a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)))*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)
     $     +cos(b))-(y3+a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))*(0.5d0*sin(b)*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*(2*y3+4*a)+sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(1+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))-0.5d0*(y1
     $     *cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     *(1+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))*(2*y3+4*a)+(y1
     $     *cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *(a/(y1**2+y2**2+(y3+2*a)**2)-a*(y3+2*a)/(y1**2+y2**2+(y3+2
     $     *a)**2)**2*(2*y3+4*a))+0.5d0/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(y2**2*cos(b)*sin(b)-a*(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*(2*y3+4*a)+1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(y2**2*cos(b)*sin(b)-a*(y1*cos(b)
     $     +(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a))*(0.5d0/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*y3+4*a)+cos(b))-1/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(-a*sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     +0.5d0*a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(1.5d0)*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)*(2*y3+4*a)-a*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*(0.5d0/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)*(2*y3+4*a)+1))))/pi/(1-nu))+0.5d0*B1*(0.125d0*(1/(y1
     $     **2+y2**2+y3**2)**(0.5d0)-1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-cos(b)*(((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)
     $     /(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))))/pi/(1
     $     -nu)+0.125d0*y2*(-1/(y1**2+y2**2+y3**2)**(1.5d0)*y2+1/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0)*y2-cos(b)*(1/(y1**2+y2**2+y3**2)
     $     *cos(b)*y2/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     -((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y2-((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3
     $     **2)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*y2
     $     -1/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y2/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*y2+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3
     $     +2*a)/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*y2))/pi/(1-nu)+1/4
     $     *((2-2*nu)*((1-2*nu)*(-1/y1/(1+y2**2/y1**2)+1/(y1*cos(b)+(y3
     $     +2*a)*sin(b))/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)+((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2
     $     *a)*sin(b))+y2**2*cos(b))+y2**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))-2*y2**2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     /(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*cos(b))/(1
     $     +y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2))*cot(b)+1/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)+y3+2*a)*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))-y2**2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     **2*(2*nu+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y2**2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)-cos(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+y2**2
     $     *cos(b)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y2**2*cos(b)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2
     $     +y2**2+(y3+2*a)**2)**(1.5d0))+(y3+a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*(2*nu/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2
     $     *a)+a/(y1**2+y2**2+(y3+2*a)**2))-y2**2*(y3+a)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*(2*nu/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)+a/(y1**2+y2**2+(y3+2*a)**2))+y2*(y3+a)/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*(-2*nu/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *y2-2*a/(y1**2+y2**2+(y3+2*a)**2)**2*y2)+(y3+a)*cos(b)/(y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1-2*nu-((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))-a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2))-y2**2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1-2*nu-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2
     $     *a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-a
     $     *(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))-y2**2*(y3+a)*cos(b)/(y1
     $     **2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))**2*(1-2*nu-((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0))-a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))
     $     +y2*(y3+a)*cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(-1
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)*y2/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *y2+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2+2*a*(y3+2*a)/(y1**2+y2
     $     **2+(y3+2*a)**2)**2*y2))/pi/(1-nu))+0.5d0*B2*(0.125d0*((-1+2
     $     *nu)*sin(b)*(1/(y1**2+y2**2+y3**2)**(0.5d0)*y2/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-1/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b)))-y1*(-1/(y1**2+y2**2+y3**2)**(1.5d0)
     $     *y2+1/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2)+(y1*cos(b)-y3
     $     *sin(b))/(y1**2+y2**2+y3**2)*cos(b)*y2/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))-(y1*cos(b)-y3*sin(b))*((y1**2
     $     +y2**2+y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)
     $     **(1.5d0)/((y1**2+y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))
     $     *y2-(y1*cos(b)-y3*sin(b))*((y1**2+y2**2+y3**2)**(0.5d0)
     $     *cos(b)-y3)/(y1**2+y2**2+y3**2)/((y1**2+y2**2+y3**2)**(0.5d0)
     $     -y1*sin(b)-y3*cos(b))**2*y2-(y1*cos(b)+(y3+2*a)*sin(b))/(y1
     $     **2+y2**2+(y3+2*a)**2)*cos(b)*y2/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+(y1*cos(b)+(y3+2*a)
     $     *sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y2+(y1*cos(b)+(y3+2
     $     *a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*y2)/pi/(1-nu)+1/4*((
     $     -2+2*nu)*(1-2*nu)*cot(b)*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)
     $     -cos(b)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))+(2-2*nu)
     $     *y1/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)+y3+2*a)**2*(2*nu+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)*y2+(2-2*nu)*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)*a/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*y2-(2
     $     -2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*(cos(b)+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *y2-(2-2*nu)*(y1*cos(b)+(y3+2*a)*sin(b))/((y1**2+y2**2+(y3+2
     $     *a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2
     $     +(y3+2*a)**2)**(1.5d0)*y2-(y3+a)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*((1-2*nu)*cot(b)-2*nu*y1/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)+y3+2*a)-a*y1/(y1**2+y2**2+(y3+2*a)**2))*y2+(y3+a)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(2*nu*y1/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)+y3+2*a)**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*y2+2*a*y1/(y1**2+y2**2+(y3+2*a)**2)**2*y2)+(y3+a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(cos(b)*sin(b)+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)*cot(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*((2-2*nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b)))+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*(sin(b)-(y3+2*a)*(y1*cos(b)+(y3+2*a)*sin(b))/(y1**2
     $     +y2**2+(y3+2*a)**2)-(y1*cos(b)+(y3+2*a)*sin(b))*((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))))/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-(y3
     $     +a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y2*cot(b)*((2-2
     $     *nu)*cos(b)-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b)))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     *cot(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*((2-2*nu)*cos(b)
     $     -((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))*y2
     $     +((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)*cot(b)
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*(-cos(b)/(y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     *cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)
     $     +(y3+2*a)*cos(b))**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*y2)-a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)*(sin(b)-(y3+2*a)*(y1
     $     *cos(b)+(y3+2*a)*sin(b))/(y1**2+y2**2+(y3+2*a)**2)-(y1*cos(b)
     $     +(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b)))*y2+a/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*(2*(y3+2*a)*(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)**2*y2-(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(y1**2+y2**2+(y3+2*a)**2)*cos(b)*y2/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))+(y1*cos(b)
     $     +(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2**2+(y3
     $     +2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*y2+(y1*cos(b)
     $     +(y3+2*a)*sin(b))*((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)
     $     +y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*y2)))/pi/(1-nu))
     $     +0.5d0*B3*(0.125d0*((2-2*nu)*(1/(y1*cos(b)-y3*sin(b))/(1+y2
     $     **2/(y1*cos(b)-y3*sin(b))**2)+((y1**2+y2**2+y3**2)**(0.5d0)
     $     *sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2*cos(b))+y2**2/(y1**2
     $     +y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)-y3*sin(b))+y2**2
     $     *cos(b))-2*y2**2*(y1**2+y2**2+y3**2)**(0.5d0)*sin(b)/(y1*(y1
     $     *cos(b)-y3*sin(b))+y2**2*cos(b))**2*cos(b))/(1+y2**2*(y1**2
     $     +y2**2+y3**2)*sin(b)**2/(y1*(y1*cos(b)-y3*sin(b))+y2**2
     $     *cos(b))**2)-1/(y1*cos(b)+(y3+2*a)*sin(b))/(1+y2**2/(y1
     $     *cos(b)+(y3+2*a)*sin(b))**2)-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))+y2**2/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1
     $     *(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))-2*y2**2*(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))**2*cos(b))/(1+y2**2*(y1**2+y2**2+(y3+2
     $     *a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))**2))+sin(b)*(((y1**2+y2**2+y3**2)**(0.5d0)*cos(b)
     $     -y3)/(y1**2+y2**2+y3**2)**(0.5d0)/((y1**2+y2**2+y3**2)
     $     **(0.5d0)-y1*sin(b)-y3*cos(b))-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b)))+y2*sin(b)*(1/(y1**2+y2**2+y3**2)*cos(b)*y2/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))-((y1**2+y2**2+y3
     $     **2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)**(1.5d0)/((y1**2
     $     +y2**2+y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))*y2-((y1**2+y2**2
     $     +y3**2)**(0.5d0)*cos(b)-y3)/(y1**2+y2**2+y3**2)/((y1**2+y2**2
     $     +y3**2)**(0.5d0)-y1*sin(b)-y3*cos(b))**2*y2-1/(y1**2+y2**2
     $     +(y3+2*a)**2)*cos(b)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*y2+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2
     $     *a)/(y1**2+y2**2+(y3+2*a)**2)/((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))**2*y2))/pi/(1-nu)+1/4
     $     *((2-2*nu)*(-1/y1/(1+y2**2/y1**2)+1/(y1*cos(b)+(y3+2*a)
     $     *sin(b))/(1+y2**2/(y1*cos(b)+(y3+2*a)*sin(b))**2)+((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)
     $     *sin(b))+y2**2*cos(b))+y2**2/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*sin(b)/(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2
     $     *cos(b))-2*y2**2*(y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*sin(b)
     $     /(y1*(y1*cos(b)+(y3+2*a)*sin(b))+y2**2*cos(b))**2*cos(b))/(1
     $     +y2**2*(y1**2+y2**2+(y3+2*a)**2)*sin(b)**2/(y1*(y1*cos(b)+(y3
     $     +2*a)*sin(b))+y2**2*cos(b))**2))+(2-2*nu)*sin(b)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)
     $     +a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))-(2-2*nu)*y2**2*sin(b)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)-(2-2*nu)*y2**2*sin(b)/((y1**2
     $     +y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*a/(y1
     $     **2+y2**2+(y3+2*a)**2)**(1.5d0)+(y3+a)*sin(b)/(y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(1+((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)
     $     **2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2))-y2**2
     $     *(y3+a)*sin(b)/(y1**2+y2**2+(y3+2*a)**2)**(1.5d0)/((y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(1+((y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2
     $     +(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)*cos(b))*(cos(b)+a
     $     /(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3+2*a)/(y1**2+y2**2
     $     +(y3+2*a)**2))-y2**2*(y3+a)*sin(b)/(y1**2+y2**2+(y3+2*a)**2)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(1+((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3
     $     +2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))+a*(y3
     $     +2*a)/(y1**2+y2**2+(y3+2*a)**2))+y2*(y3+a)*sin(b)/(y1**2+y2
     $     **2+(y3+2*a)**2)**(0.5d0)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*(1/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)*y2/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1
     $     *sin(b)+(y3+2*a)*cos(b))*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0))-((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)*cos(b)+y3+2*a)
     $     /((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)-y1*sin(b)+(y3+2*a)
     $     *cos(b))**2*(cos(b)+a/(y1**2+y2**2+(y3+2*a)**2)**(0.5d0))/(y1
     $     **2+y2**2+(y3+2*a)**2)**(0.5d0)*y2-((y1**2+y2**2+(y3+2*a)**2)
     $     **(0.5d0)*cos(b)+y3+2*a)/((y1**2+y2**2+(y3+2*a)**2)**(0.5d0)
     $     -y1*sin(b)+(y3+2*a)*cos(b))*a/(y1**2+y2**2+(y3+2*a)**2)
     $     **(1.5d0)*y2-2*a*(y3+2*a)/(y1**2+y2**2+(y3+2*a)**2)**2*y2))
     $     /pi/(1-nu))

        return
        end
        
