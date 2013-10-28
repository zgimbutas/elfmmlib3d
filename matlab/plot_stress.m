figure(11)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(1,1,:))
axis equal
title('\sigma_{13}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(21)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(2,1,:))
axis equal
title('\sigma_{23}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(31)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(3,1,:))
axis equal
title('\sigma_{33}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar


figure(12)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(1,2,:))
axis equal
title('\sigma_{12}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(22)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(2,2,:))
axis equal
title('\sigma_{22}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(32)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(3,2,:))
axis equal
title('\sigma_{32}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar



figure(13)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(1,3,:))
axis equal
title('\sigma_{13}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(23)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(2,3,:))
axis equal
title('\sigma_{23}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(33)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.stress(3,3,:))
axis equal
title('\sigma_{33}')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar


figure(4)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', ...
   (U.stress(1,1,:)+U.stress(2,2,:)+U.stress(3,3,:))/3)
axis equal
title('\sigma_{ii}/3')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar