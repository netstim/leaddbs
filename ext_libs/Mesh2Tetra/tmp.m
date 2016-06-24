 compile_c_files

 % Load the Training Data
 load('example_jaw');

 % Do normal delaunay tetrahedron generation
 Tn = delaunayn(FV.vertices);

 % Do constrained tetrahedron generation
 Tc = Mesh2Tetra(FV.vertices,FV.faces);

 figure,
 subplot(2,2,1), hold on;
  patch(FV,'facecolor',[1 1 0]);
  axis equal; title('Surface Mesh');
  plot3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),'r*');
 subplot(2,2,2), hold on;
  tetramesh(Tn,FV.vertices,'facecolor',[0 1 0]);
  axis equal; title('Tetrahedrons After Delaunayn');
 subplot(2,2,3), hold on;
  tetramesh(Tc,FV.vertices,'facecolor',[0 0 1]);
  axis equal; title('Tetrahedrons After Mesh2Tetra');