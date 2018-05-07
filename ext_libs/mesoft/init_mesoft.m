% mex all the ccode
cd(fullfile(ea_getearoot,'ext_libs','mesoft','ccode'));
%display('mexing c-code');
%ft_make
cd ..

% add the neccesary paths
disp('adding MATLAB-path');
addpath(pwd)
addpath(fullfile(pwd,'ccode'));
addpath(fullfile(pwd,'aux_mfiles'));
addpath(fullfile(pwd,'sphereInterpolation'));