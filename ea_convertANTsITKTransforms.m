function ea_convertANTsITKTransforms(inputFile, outputFile, invert)
% ea_convertANTsITKTransform - read an ITK (.txt, .mat, .tfm) transform file
% and and convert it to .mat format using antsApplyTransforms. Invert if
% needed.
%
% No further conversion are applied here! To be useful in matlab a conversion of
% the matrix to RAS might be needed, depending on the application.
%
% Andreas Husch, University of Luxembourg, Interventional Neuroscience
% Group, 2019

if(nargin < 3)
    invert = 0;
end

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransforms = ea_path_helper([basedir, 'ext_libs', filesep, 'ANTs', filesep, 'antsApplyTransforms.exe']);
else
    applyTransforms = [basedir, 'ext_libs', filesep, 'ANTs', filesep, 'antsApplyTransforms.', computer('arch')];
end

cmd = [applyTransforms ...
    ' -t ' inputFile ' -o Linear[' outputFile ', ' num2str(invert) ']' ];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

end