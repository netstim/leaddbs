if ismac
   compflags = '';
elseif isunix
   compflags = 'COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
elseif ispc
   compflags = 'COMPFLAGS="$COMPFLAGS /MT"';
end

mex(compflags, 'AccumulateBilin.cpp');
mex(compflags, 'AccumulateBilinWeighted.cpp');
mex(compflags, 'anisoDiffusion.cpp');
mex(compflags, 'anisoDiffusionHomogenous.cpp');
mex('-compatibleArrayDims', compflags, 'BuildFibres.cpp');
mex('-compatibleArrayDims', compflags, 'CreateConnectivityMatrixROI.cpp');
mex('-compatibleArrayDims', compflags, 'pcRJMCMC.cpp');
mex('-compatibleArrayDims', compflags, 'reparametrize_arclen.cpp');
