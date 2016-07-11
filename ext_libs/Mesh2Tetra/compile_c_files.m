% This script will compile all the C files of the registration methods
cd(['functions', filesep, 'mexfunctions']);

mex BarycentricCoordinatesTetrahedron.c -v                         
mex BarycentricCoordinatesTriangle.c -v                             
mex CheckInsideFace.c -v                            
mex LineLineIntersect.c -v              
mex LineTriangleIntersection.c -v  
mex SphereFrom4Points.c -v                        
mex TriangleTriangleIntersection.c -v
mex CheckVolumeFaceMesh.c -v
mex CheckVolumeTetraMesh.c -v

cd(['..', filesep, '..']);
