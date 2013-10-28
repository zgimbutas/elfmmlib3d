function [U]=elh3dtriadirect(rlam,rmu,nsource,triaflat,trianorm,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifstrain,ntarget,target,ifpottarg,ifstraintarg)
%ELHFMM3DTRIADIRECT Elastostatic triangle interactions in R^3, direct algorithm.
%
% Half space Green's function with zero normal stress at z=0
%
%
% [U]=ELH3DTRIADIRECT(RLAM,RMU,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,IFPOT,IFSTRAIN);
%
% [U]=ELH3DTRIADIRECT(RLAM,RMU,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,IFPOT,IFSTRAIN,...
%         NTARGET,TARGET,IFPOTTARG,IFSTRAINTARG);
%
%
% This subroutine evaluates the elastostatic potential and strain due
% to a collection of flat triangles with constant single and/or
% double layer densities. We use 
% ((2-A) delta_ij 1/r+ A r_i r_j/r^3)/(2*rmu), A=(rlam+rmu)/(rlam+2*rmu)
% for the Green's function, without the (1/4 pi) scaling.  
% Self-interactions are included.
%
% It is capable of evaluating the layer potentials either on or 
% off the surface (or both).            
%
%
% Input parameters:
% 
% rlam, rmu - Lame parameters
% nsource - number of triangles
% triaflat - double (3,3,ntriangles): array of triangle vertex coordinates
% trianorm - double (3,ntriangles): triangle normals
% source - double (3,ntriangles): triangle centroids
% ifcharge - single layer computation flag
%
%         0 => do not compute
%         1 => include elastostatic SLP contribution
% 
% charge - complex (3,ntriangles): piecewise constant SLP (charge) strength 
% ifdipole - double layer computation flag
%
%         0 => do not compute
%         1 => include elastostatic DLP contribution
% 
% dipole - complex (3,ntriangles): piecewise constant DLP (dipole) strength 
%
%     In the present version, dipole orientation vector is assumed to
%     BE SET EQUAL to the triangle normal. 
%
% ifpot - displacement computation flag, 
%         1 => compute the displacement, otherwise no
% ifstrain - strain computation flag, 1 => compute the strain, otherwise no
%
% ntarget - number of targets
% target - double (3,ntarget): target locations
%
% ifpottarg - target displacement computation flag, 
%         1 => compute the displacement, otherwise no
% ifstraintarg - target strain computation flag, 
%         1 => compute the strain, otherwise no
%
%
% Output parameters: 
%
% U.pot - complex (3nsource) - displacement at triangle centroids
% U.strain - complex (3,3,nsource) - strain at triangle centroids
% U.pottarg - complex (3,ntarget) - displacement at targets
% U.straintarg - complex (3,3,ntarget) - strain at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%

if( nargin == 10 ) 
  ifpot = 1;
  ifstrain = 1;
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifstraintarg = 0;
end

if( nargin == 12 ) 
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifstraintarg = 0;
end

if( nargin == 14 ) 
  ifpottarg = 1;
  ifstraintarg = 1;
end

pot=zeros(3,nsource);
strain=zeros(3,3,nsource);
pottarg=zeros(3,ntarget);
straintarg=zeros(3,3,ntarget);
ier=0;


mex_id_ = 'elh3dtriadirecttarg(i double[x], i double[x], i double[], i double[xx], i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i int[x], io double[], i int[x], io double[], i int[x], i double[xx], i int[x], io double[], i int[x], io double[])';
[pot, strain, pottarg, straintarg] = elfmm3d_r2013a(mex_id_, rlam, rmu, triaflat, trianorm, nsource, source, ifcharge, charge, ifdipole, dipstr, ifpot, pot, ifstrain, strain, ntarget, target, ifpottarg, pottarg, ifstraintarg, straintarg, 1, 1, 3, nsource, 1, 3, nsource, 1, 1, 1, 1, 1, 3, ntarget, 1, 1);

if( ifpot == 1 ) U.pot=pot; end
if( ifstrain == 1 ) U.strain=reshape(strain,3,3,nsource); end
if( ifpottarg == 1 ) U.pottarg=pottarg; end
if( ifstraintarg == 1 ) U.straintarg=reshape(straintarg,3,3,ntarget); end
U.ier=ier;




