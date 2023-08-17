if ismac
    compflags = '';
elseif isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
    % compflags = [' COMPFLAGS=''$COMPFLAGS -static-libstdc++''', ...
    %             ' -DPARALLEL_OPENMP', ...
    %             ' CXXFLAGS=''$CXXFLAGS -fopenmp''', ...
    %             ' LINKLIBS=''$LINKLIBS -lgomp'''];
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT"';
end

mex(compflags, 'AccumulateTrilin.cpp');
mex(compflags, 'AccumulateTrilinWeighted.cpp');
mex(compflags, 'BuildFibres.cpp');
mex(compflags, 'pcRJMCMC.cpp');
mex(compflags, 'printTOstderr.cpp');
mex(compflags, 'reparametrize_arclen.cpp');
mex(compflags, 'SelectCorticalFibers.cpp');
