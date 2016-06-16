% mex all the ccode
%display('mexing c-code');
cd ccode
%ft_make
cd ..

% add the neccesary paths
display('adding MATLAB-path');
addpath(pwd)
addpath(fullfile(pwd,'ccode'));
addpath(fullfile(pwd,'aux_mfiles'));
addpath(fullfile(pwd,'sphereInterpolation'));