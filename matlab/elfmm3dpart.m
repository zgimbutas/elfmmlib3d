function [U]=elfmm3dpart(iprec,rlam,rmu,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifstrain,ntarget,target,ifpottarg,ifstraintarg)
%ELFMM3DPART Elastostatic particle target FMM in R^3.
%
% Elastostatic FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interaction) and interactions with targets.
%
% [U]=ELFMM3DPART(IPREC,RLAM,RMU,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFSTRAIN);
%
% [U]=ELFMM3DPART(IPREC,RLAM,RMU,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFSTRAIN,...
%         NTARGET,TARGET,IFPOTTARG,IFSTRAINTARG);
%
%
% This subroutine evaluates the elastostatic potential and strain due
% to a collection of elastostatic charges and dipoles. 
% ((2-A) delta_ij 1/r+ A r_i r_j/r^3)/(2*rmu), A=(rlam+rmu)/(rlam+2*rmu)
% for the Green's function, without the (1/4 pi) scaling.  
% Self-interactions are not-included.
%
%
% Input parameters:
% 
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
%
% rlam, rmu - Lame parameters
% nsource - number of sources
% source - double (3,nsource): source locations
% ifcharge - charge computation flag
%
%         0 => do not compute
%         1 => include elastostatic charge contribution
% 
% charge - complex (3,nsource): charge strengths
% ifdipole - dipole computation flag
%
%         0 => do not compute
%         1 => include elastostatic dipole contribution
% 
% dipole - complex (3,nsource): dipole strengths
% dipvec - double (3,nsource): dipole orientation vectors 
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
% U.pot - complex (3nsource) - displacement at source locations
% U.strain - complex (3,3,nsource) - strain at source locations
% U.pottarg - complex (3,ntarget) - displacement at targets
% U.straintarg - complex (3,3,ntarget) - strain at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
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
if( ntarget > 0 ),
pottarg=zeros(3,ntarget);
straintarg=zeros(3,3,ntarget);
end
ier=0;

if( ntarget == 0 ) 
mex_id_ = 'elfmm3dpart(io int[x], i int[x], i double[x], i double[x], i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i double[xx], i int[x], io double[], i int[x], io double[])';
[ier, pot, strain] = elfmm3d_r2013a(mex_id_, ier, iprec, rlam, rmu, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, ifstrain, strain, 1, 1, 1, 1, 1, 3, nsource, 1, 1, 3, nsource, 1, 1);
else
mex_id_ = 'elfmm3dparttarg(io int[x], i int[x], i double[x], i double[x], i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i double[xx], i int[x], io double[], i int[x], io double[], i int[x], i double[xx], i int[x], io double[], i int[x], io double[])';
[ier, pot, strain, pottarg, straintarg] = elfmm3d_r2013a(mex_id_, ier, iprec, rlam, rmu, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, ifstrain, strain, ntarget, target, ifpottarg, pottarg, ifstraintarg, straintarg, 1, 1, 1, 1, 1, 3, nsource, 1, 1, 3, nsource, 1, 1, 1, 3, ntarget, 1, 1);
end

if( ifpot == 1 ) U.pot=pot; end
if( ifstrain == 1 ) U.strain=reshape(strain,3,3,nsource); end
if( ifpottarg == 1 ) U.pottarg=pottarg; end
if( ifstraintarg == 1 ) U.straintarg=reshape(straintarg,3,3,ntarget); end
U.ier=ier;


