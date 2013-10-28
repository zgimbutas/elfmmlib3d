figure(41)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.pot(1,:))
axis equal
title('d_1')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(42)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.pot(2,:))
axis equal
title('d_2')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(43)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', U.pot(3,:))
axis equal
title('d_3')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar