% mex all the ccode
%display('mexing c-code');
cd(fullfile(ea_getearoot,'ext_libs','mesoft','ccode'));
%ft_make
cd ..

% add the neccesary paths
display('adding MATLAB-path');
addpath(pwd)
addpath(fullfile(pwd,'ccode'));
addpath(fullfile(pwd,'aux_mfiles'));
addpath(fullfile(pwd,'sphereInterpolation'));