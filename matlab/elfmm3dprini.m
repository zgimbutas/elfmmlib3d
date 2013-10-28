function elfmm3dprini(unit1,unit2)
%PRINI Initialize simple printing routines.
%
% Calling PRINI(6,13) causes printing to screen and file fort.13.     
%

if (nargin == 1 )
unit2=0;
end

mex_id_ = 'prini(i int[x], i int[x])';
elfmm3d_r2013a(mex_id_, unit1, unit2, 1, 1);

