% Set up a random triangle
x = -rand(1,3);
y = -rand(1,3);
z = -rand(1,3);

ns = 1000
sx=-rand(ns,1);
sy=-rand(ns,1);
sz=-rand(ns,1);

pr = .45;

ss = .2;
ts = .1;
ds = .3;

'Displacements'
tic
U=CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
total_time = toc
speed = ns/total_time

'Strains'
tic
S=CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
total_time = toc
speed = ns/total_time

