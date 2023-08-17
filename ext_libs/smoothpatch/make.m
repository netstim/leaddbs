if ismac
    compflags = '';
elseif isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT"';
end

mex(compflags, 'ea_smoothpatch_curvature_double.c');
mex(compflags, 'ea_smoothpatch_inversedistance_double.c');
mex(compflags, 'ea_vertex_neighbours_double.c');