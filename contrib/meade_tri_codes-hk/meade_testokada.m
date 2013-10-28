%
%  Material constants
%

pr = .3
rmu = 2.3
rlam = 2*rmu*pr/(1-2*pr)

%pr = .25
%rmu = 1
%rlam = 2*rmu*pr/(1-2*pr)


PAGER('cat')
format long

sx = .2 
sy = .3
sz = 0

ss = 0;
ds = 0;
ts = 1;


% Construct a unit triangle
x = [0 1 1]
y = [0 0 1]
z = -[-1 -1 -1]

Displacement = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
Strain = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
Stress = StrainToStress(Strain, rlam, rmu);

u1 = [Displacement.x Displacement.y Displacement.z];
s1 = [Strain.xx Strain.xy Strain.xz;
      Strain.xy Strain.yy Strain.yz;
      Strain.xz Strain.yz Strain.zz];

% Construct a unit triangle
x = [1 0 0]
y = [1 1 0]
z = -[-1 -1 -1]

Displacement = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
Strain = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
Stress = StrainToStress(Strain, rlam, rmu);

u2 = [Displacement.x Displacement.y Displacement.z];
s2 = [Strain.xx Strain.xy Strain.xz;
      Strain.xy Strain.yy Strain.yz;
      Strain.xz Strain.yz Strain.zz];


u1+u2
s1+s2
