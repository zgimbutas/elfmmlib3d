figure(1)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', dipstr(1,:))
axis equal
title('u_1')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(2)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', dipstr(2,:))
axis equal
title('u_2')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar
figure(3)
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)', dipstr(3,:))
axis equal
title('u_3')
xlabel('x'); ylabel('y'); zlabel('z');
colorbar