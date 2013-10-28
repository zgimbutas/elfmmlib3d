%
% Triangulation of a unit square
%

verts = zeros(3,4);
iverts = zeros(3,2);

verts(1:3,1) = [0;0;0];
verts(1:3,2) = [1;0;0];
verts(1:3,3) = [1;1;0];
verts(1:3,4) = [0;1;0];

iverts(1:3,1) = [1;2;3];
iverts(1:3,2) = [1;3;4];

[ndim,nverts] = size(verts);

if( 2 == 2 )
%
% shift and duplicate
%
scale = [1;1;1];
shift = [1;0;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [0;1;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [2;0;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [0;2;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [4;0;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [0;4;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [8;0;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [0;8;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1];

[ndim,nverts] = size(verts);

%
% shift and duplicate
%
scale = [1;1;1];
shift = [16;0;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);

end

%%%verts(3,:) =  sin(1+.2*verts(1,:)) + sin(.3*verts(2,:));
verts(3,:) =  sin(1+.2*verts(1,:));
%%%verts(3,:) =  sin(1+.1*verts(1,:));


verts0=verts;
iverts0=iverts;
%
% shift and duplicate
%
scale_vec = [.5 .5 .5 .5  .5
              2  1  1  1   1.5
              1  1  1  1   1];
shift_vec = [10 30  15 3    5
             20 60  80 -20  50
             -1  1  -1 0    4];

n=5;

for i=1:n
scale = scale_vec(1:3,i);
shift = shift_vec(1:3,i);

verts1 = [verts0.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts0 + i*nverts;

verts = [verts verts1];
iverts = [iverts iverts1];
end

perm=[2  3  4  1  5];

for i=1:n
scale = scale_vec(1:3,i);
shift = shift_vec(1:3,perm(i));
shift = shift + [10; -10; -10];

verts1 = [verts0.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts0 + (n+i)*nverts;

verts = [verts verts1];
iverts = [iverts iverts1];
end


[ndim,nverts] = size(verts);

if( 1 == 2 )
%
% shift and duplicate
%
scale = [1;1;1];
%%shift = [0;32;0];
shift = [0;120;0];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);
end

if( 1 == 2 )
%
% shift and duplicate
%
scale = [1;1;1];
%%shift = [0;0;32];
shift = [0;0;60];

verts1 = [verts.*repmat(scale,1,nverts)+repmat(shift,1,nverts)];
iverts1 = iverts + nverts;

verts = [verts verts1];
iverts = [iverts iverts1]

[ndim,nverts] = size(verts);
end


%
% scale
%
verts = verts/16;


%
%  rotate the geometry around the z-axis by angle alpha
%
alpha = 0;
rotmat = [cos(alpha) -sin(alpha) 0; 
          sin(alpha)  cos(alpha) 0;
             0             0     1];
verts = rotmat * verts;

%
%  rotate the geometry around the y-axis by angle beta
%
beta = -pi/2 + 1e-6;
rotmat = [cos(beta) 0 -sin(beta); 
          0         1          0;
          sin(beta) 0  cos(beta)];
verts = rotmat * verts;

%figure(1)
%plot3(verts(1,:),verts(2,:),verts(3,:),'*')
%axis equal

if( 2 == 2 )
[ndim,nverts] = size(verts);
[ndim,nfaces] = size(iverts);
ifaces = iverts;

figure(2)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)')
axis equal
end