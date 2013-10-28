%
%  Test elastostatic triangle FMMs in R^3
%

%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere180.a.tri';
%filename_geo = 'sphere320.a.tri';
%filename_geo = 'sphere1280.a.tri';
%filename_geo = 'sphere2880.a.tri';
%filename_geo = 'sphere11520.a.tri';
%filename_geo = 'sphere20480.a.tri';

fid = fopen(filename_geo,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[3,nfaces]);

fclose(fid);


%
% For half space problem, z component must be less or equal to zero
%
verts(3,:)=verts(3,:) - 5;

%
%  refined rectangle
%
%test2a
%filename_geo='rectangle'

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

nsource = ntri

%
%  timings
%

ifcharge=1;
charge = ones(3,ntri);
ifdipole=1;
dipstr = ones(3,ntri);
dipvec = trianorm;

ifcharge=1;
charge = rand(3,ntri);
ifdipole=1;
dipstr = rand(3,ntri);
dipvec = trianorm;

%%dipstr = zeros(3,ntri);
%%dipstr(2,:) = 1;
%%dipstr(2,:) = 1*cos(.1*source(3,:));
%%dipstr = cross(dipstr,cross(dipstr,trianorm));

printf('======================================\n')
ifcharge
ifdipole
ifpot = 1
ifstrain = 1

ntarget = ntri
target = source(:,1:ntri);
target(1,:) = target(1,:) + 10;

%%h=1e-4
%%target = [source + h*trianorm source - h*trianorm];
%%target = [source + h*trianorm];
%%target = [source - h*trianorm];

[ndim,ntarget] = size(target);
ifpottarg = 1
ifstraintarg = 1
ntarget
%ntarget = 0


'Lame parameters'

rlam=3.2
rmu=2.3


'Elastostatic triangle target FMM in R^3 (half space)'

tic
iprec=1
[U]=elhfmm3dtria(iprec,rlam,rmu,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifstrain,ntarget,target,ifpottarg,ifstraintarg);
total_time=toc
speed=(ntri+ntarget)/total_time

'Elastostatic triangle target direct evaluation in R^3 (half space)'
tic
[F]=elh3dtriadirect(rlam,rmu,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifstrain,ntarget,target,ifpottarg,ifstraintarg);
total_time=toc
speed=(ntri+ntarget)/total_time


if( ifpot ), U.pot=U.pot/(4*pi); end
if( ifstrain ), U.strain=U.strain/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( ifstrain ), F.strain=F.strain/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
end

if( ifstrain ),
%rms_strain = norm(reshape(F.strain,9,nsource),2)/sqrt(nsource)
rms_error_strain = norm(reshape(U.strain - F.strain,9,nsource),2)/sqrt(nsource)
end
%%%break;

if( ifpottarg ), U.pottarg=U.pottarg/(4*pi); end
if( ifstraintarg ), U.straintarg=U.straintarg/(4*pi); end

if( ifpottarg ), F.pottarg=F.pottarg/(4*pi); end
if( ifstraintarg ), F.straintarg=F.straintarg/(4*pi); end

if( ifpottarg ),
%rms_pottarg = norm((F.pottarg),2)/sqrt(nsource)
rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(ntarget)
end

if( ifstraintarg ),
rms_straintarg = norm(reshape(F.straintarg,9,ntarget),2)/sqrt(ntarget)
rms_error_straintarg = ...
    norm(reshape(U.straintarg - F.straintarg,9,ntarget),2)/sqrt(ntarget)
end
%%%break;



'Elastostatic triangle FMM in R^3 (half space)'

ifpottarg = 0
ifstraintarg = 0


tic
iprec=1
[U]=elhfmm3dtria(iprec,rlam,rmu,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifstrain);
total_time=toc
speed=(ntri)/total_time

'Elastostatic triangle direct evaluation in R^3 (half space)'

tic
[F]=elh3dtriadirect(rlam,rmu,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifstrain);
total_time=toc
speed=(ntri)/total_time


if( ifpot ), U.pot=U.pot/(4*pi); end
if( ifstrain ), U.strain=U.strain/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( ifstrain ), F.strain=F.strain/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
end

if( ifstrain ),
%rms_strain = norm(reshape(F.strain,9,nsource),2)/sqrt(nsource)
rms_error_strain = norm(reshape(U.strain - F.strain,9,nsource),2)/sqrt(nsource)
end
%%%break;


%filename_out=['output_' filename_geo '.mat']
%save(filename_out,'-v6')

%plot_solution
%plot_displacement
%plot_slip


