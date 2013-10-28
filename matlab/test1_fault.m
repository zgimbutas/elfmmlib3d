%
%  Test elastostatic triangle FMMs in R^3
%

%
%  Refined faults
%
test1_geom
filename_geo='faults11demo'

nverts,nfaces

%
%  create triangle vertex and normal arrays
%

ntri = nfaces;
triangles = zeros(3,3,ntri);

for i=1:ntri
	
%triangles(1:3,1,i) = verts(1:3,ifaces(1,i));
%triangles(1:3,2,i) = verts(1:3,ifaces(2,i));
%triangles(1:3,3,i) = verts(1:3,ifaces(3,i));

triangles(1:3,1:3,i) = verts(1:3,ifaces(1:3,i));

end

%
%  build triangle normals and area vector
%
trianorm = zeros(3,ntri);
triaarea = zeros(1,ntri);

for i=1:ntri

vec1 = triangles(1:3,2,i) - triangles(1:3,1,i);
vec2 = triangles(1:3,3,i) - triangles(1:3,1,i);

trianorm(1:3,i) = cross(vec1,vec2);
triaarea(i) = norm(trianorm(1:3,i))/2;

trianorm(1:3,i) = trianorm(1:3,i)/norm(trianorm(1:3,i));

end

%%% sum(triaarea), 4*pi

%
%  centroids
%

source = sum(triangles,2)/3;
source = reshape(source,3,ntri);

%plot3(source(1,:),source(2,:),source(3,:))

%
%  timings
%

ifcharge=0;
charge = ones(3,ntri);
ifdipole=1;
dipstr = ones(3,ntri);
dipvec = trianorm;

dipstr = zeros(3,ntri);
dipstr(2,:) = 1;
%%dipstr(2,:) = 1*cos(.1*source(3,:));
%%dipstr = cross(dipstr,cross(dipstr,trianorm));



ifcharge
ifdipole
ifpot = 1
ifstrain = 1

ntarget = min(10,ntri)
target = source(:,ntri);
target(1,:) = target(1,:) + 3;

h=1e-4
target = [source + h*trianorm source - h*trianorm];
target = [source + h*trianorm];
target = [source - h*trianorm];
[ndim,ntarget] = size(target)
ifpottarg = 0
ifstraintarg = 0
ntarget = 0
'Elastostatic triangle target FMM in R^3'

rlam=3.2
rmu=2.3

tic
iprec=1;
[U]=elfmm3dtria(iprec,rlam,rmu,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifstrain,ntarget,target,ifpottarg,ifstraintarg);
total_time=toc

U.pot=U.pot/(4*pi);
U.strain=U.strain/(4*pi);


U.stress = 2*rmu*U.strain;
d = (U.strain(1,1,:)+U.strain(2,2,:)+U.strain(3,3,:));
U.stress(1,1,:) = U.stress(1,1,:) + rlam*d;
U.stress(2,2,:) = U.stress(2,2,:) + rlam*d;
U.stress(3,3,:) = U.stress(3,3,:) + rlam*d;


%filename_out=['output_' filename_geo '.mat']
%save(filename_out,'-v6')


plot_slip

plot_strain
%%plot_stress

plot_displacement

