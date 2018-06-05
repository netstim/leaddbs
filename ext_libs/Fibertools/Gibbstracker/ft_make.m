if ismac
    compflags = '';
%     compflags = ' CXXFLAGS=''$CXXFLAGS -Wno-deprecated-register''';
elseif isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT"';
end

eval(['mex' compflags, ' AccumulateBilin.cpp']);
eval(['mex' compflags, ' AccumulateBilinWeighted.cpp']);
eval(['mex' compflags, ' anisoDiffusion.cpp']);
eval(['mex' compflags, ' anisoDiffusionHomogenous.cpp']);
eval(['mex' compflags, ' BuildFibres.cpp']);
eval(['mex' compflags, ' CreateConnectivityMatrixROI.cpp']);
eval(['mex' compflags, ' pcRJMCMC.cpp']);
eval(['mex' compflags, ' reparametrize_arclen.cpp']);
