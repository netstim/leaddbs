% This script will compile all the C files of the registration methods
cd(['functions', filesep, 'mexfunctions']);

if ismac
    compflags = '';
elseif isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT"';
end

mex('-compatibleArrayDims', compflags, 'BarycentricCoordinatesTetrahedron.c');
mex('-compatibleArrayDims', compflags, 'BarycentricCoordinatesTriangle.c');                         
mex('-compatibleArrayDims', compflags, 'CheckInsideFace.c');
mex('-compatibleArrayDims', compflags, 'LineLineIntersect.c');
mex('-compatibleArrayDims', compflags, 'LineTriangleIntersection.c');
mex('-compatibleArrayDims', compflags, 'SphereFrom4Points.c');
mex('-compatibleArrayDims', compflags, 'TriangleTriangleIntersection.c');
mex('-compatibleArrayDims', compflags, 'CheckVolumeFaceMesh.c');
mex('-compatibleArrayDims', compflags, 'CheckVolumeTetraMesh.c');

cd(['..', filesep, '..']);
