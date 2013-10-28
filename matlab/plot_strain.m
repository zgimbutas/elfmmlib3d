figure(11)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(1,1,:))
axis equal
title('\epsilon_{13}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(21)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(2,1,:))
axis equal
title('\epsilon_{23}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(31)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(3,1,:))
axis equal
title('\epsilon_{33}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar


figure(12)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(1,2,:))
axis equal
title('\epsilon_{12}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(22)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(2,2,:))
axis equal
title('\epsilon_{22}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(32)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(3,2,:))
axis equal
title('\epsilon_{32}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar



figure(13)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(1,3,:))
axis equal
title('\epsilon_{13}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(23)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(2,3,:))
axis equal
title('\epsilon_{23}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(33)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.strain(3,3,:))
axis equal
title('\epsilon_{33}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar


figure(4)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', ...
   (U.strain(1,1,:)+U.strain(2,2,:)+U.strain(3,3,:))/3)
axis equal
title('\epsilon_{ii}/3')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar