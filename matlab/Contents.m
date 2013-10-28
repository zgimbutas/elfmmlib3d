% Elastostatic FMMs in R^3.
%
% Triangle FMM routines (constant densities on flat triangles).
%   elfmm3dtria  - Elastostatic triangle FMM in R^3.
%   elhfmm3dtria  - Elastostatic triangle FMM in half space.
%
% Particle FMM routines.
%   elfmm3dpart  - Elastostatic particle FMM in R^3.
%   elhfmm3dpart  - Elastostatic particle FMM in half space.
%
% Direct evaluation routines (constant densities on flat triangles).
%   el3dtriadirect  - Elastostatic triangle interactions in R^3.
%   elh3dtriadirect  - Elastostatic triangle interactions in half space.
%
% Direct evaluation routines (particles).
%   el3dpartdirect  - Elastostatic particle interactions in R^3.
%   elh3dpartdirect  - Elastostatic particle interactions in half space.
%
% Triangulations.
%   atriread - Retrieve Cart3d triangulation from a file. (flat)
%   atriwrite - Store Cart3d triangulation to a file. (flat)
%   atriproc - Process triangulations in Cart3d format. (flat)
%   atrirefine - Refine Cart3d triangulation. (flat)
%   atriplot - Plot Cart3d triangulation. (flat)
%
% Utility functions. 
%   el3dstrain2stress - convert strain tensor to stress tensor.
%   el3dstrain2traction - convert strain tensor to traction vector.
%   el3dstress2traction - convert stress tensor to traction vector.
%
% Iterative methods.
%   gmres_simple    - GMRES algorithm.
%   bicgstab_simple - BiCG(stab) algorithm.
%
% Demos.
%   test_el3dpart_direct - test elastostatic particle FMM vs direct.
%   test_el3dtria_direct - test elastostatic triangle FMM vs direct.
%   test_elh3dpart_direct - test elastostatic half.sp. particle FMM vs direct.
%   test_elh3dtria_direct - test elastostatic half.sp. triangle FMM vs direct.
%   test1_fault - test elastostatic triangle FMM for faults in free space.
%
% Internal utility functions.
%   elfmm3dprini   - initialize simple printing routines.
%

%% Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
%% Contact: greengard@cims.nyu.edu
%% 
%% This program is free software; you can redistribute it and/or modify 
%% it under the terms of the GNU General Public License as published by 
%% the Free Software Foundation; either version 2 of the License, or 
%% (at your option) any later version.  This program is distributed in 
%% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
%% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%% PARTICULAR PURPOSE.  See the GNU General Public License for more 
%% details. You should have received a copy of the GNU General Public 
%% License along with this program; 
%% if not, see <http://www.gnu.org/licenses/>.
