if ismac
    compflags = '';
%     compflags = ' CXXFLAGS=''$CXXFLAGS -Wno-deprecated-register''';
elseif isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
%     compflags = [' COMPFLAGS=''$COMPFLAGS -static-libstdc++''', ...
%                 ' -DPARALLEL_OPENMP', ...
%                 ' CXXFLAGS=''$CXXFLAGS -fopenmp''', ...
%                 ' LINKLIBS=''$LINKLIBS -lgomp'''];
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT"';
end

eval(['mex' compflags, ' AccumulateTrilin.cpp']);
eval(['mex' compflags, ' AccumulateTrilinWeighted.cpp']);
eval(['mex' compflags, ' BuildFibres.cpp']);
eval(['mex' compflags, ' pcRJMCMC.cpp']);
eval(['mex' compflags, ' printTOstderr.cpp']);
eval(['mex' compflags, ' reparametrize_arclen.cpp']);
eval(['mex' compflags, ' SelectCorticalFibers.cpp']);
