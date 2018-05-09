if ismac
    compflgs = '';
elseif isunix
    compflgs = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++''';
elseif ispc
    compflgs = ' COMPFLAGS="$COMPFLAGS /MT"';
end

eval(['mex' compflgs, ' AccumulateTrilin.cpp']);
eval(['mex' compflgs, ' AccumulateTrilinWeighted.cpp']);
eval(['mex' compflgs, ' BuildFibres.cpp']);
eval(['mex' compflgs, ' pcRJMCMC.cpp']);
eval(['mex' compflgs, ' printTOstderr.cpp']);
eval(['mex' compflgs, ' reparametrize_arclen.cpp']);
eval(['mex' compflgs, ' SelectCorticalFibers.cpp']);
