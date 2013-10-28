% Construct a random triangle
x = -rand(1,3);
y = -rand(1,3);
z = -rand(1,3);

sx=-rand(4,1);
sy=-rand(4,1);
sz=-rand(4,1);

pr = .45;

ss = .2;
ts = .1;
ds = .3;

U=CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
Strain=CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);

