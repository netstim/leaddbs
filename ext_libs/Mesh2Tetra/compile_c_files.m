% This script will compile all the C files of the registration methods
cd(['functions', filesep, 'mexfunctions']);

if ismac
    compflags = '';
elseif isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT"';
end

eval(['mex' compflags, ' BarycentricCoordinatesTetrahedron.c']);
eval(['mex' compflags, ' BarycentricCoordinatesTriangle.c']);                         
eval(['mex' compflags, ' CheckInsideFace.c']);
eval(['mex' compflags, ' LineLineIntersect.c']);
eval(['mex' compflags, ' LineTriangleIntersection.c']);
eval(['mex' compflags, ' SphereFrom4Points.c']);
eval(['mex' compflags, ' TriangleTriangleIntersection.c']);
eval(['mex' compflags, ' CheckVolumeFaceMesh.c']);
eval(['mex' compflags, ' CheckVolumeTetraMesh.c']);

cd(['..', filesep, '..']);
