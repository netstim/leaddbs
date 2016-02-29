function ea_brainsfit(varargin)
% Wrapper for BRAINSFit

fixedVolume=varargin{1};
movingVolume=varargin{2};
outputVolume=varargin{3};
if nargin>3
    writematout=varargin{4};
else
    writematout=1;
end
if fileparts(movingVolume)
    volumedir = [fileparts(movingVolume), filesep];
else
    volumedir =['.', filesep];
end

fixedVolume = ea_path_helper(fixedVolume);
movingVolume = ea_path_helper(movingVolume);
outputVolume = ea_path_helper(outputVolume);

fixparams = [' --fixedVolume ' , fixedVolume, ...
             ' --movingVolume ', movingVolume, ...
             ' --outputVolume ', outputVolume, ...
             ' --useRigid --useAffine' ...
             ' --samplingPercentage 0.005' ...
             ' --removeIntensityOutliers 0.005' ...
             ' --interpolationMode Linear' ...
             ' --outputTransform ', ea_path_helper([volumedir, 'ct2anat.txt'])];

% first attempt...
paramset{1} = [fixparams, ...
               ' --initializeTransformMode useGeometryAlign' ...
               ' --maskProcessingMode ROIAUTO' ...
               ' --ROIAutoDilateSize 3'];
% second attempt...
paramset{2} = [fixparams, ...
               ' --initializeTransformMode useGeometryAlign'];          
% third attempt...
if exist([volumedir, 'ct2anat.txt'],'file')
    paramset{3} = [fixparams, ...
                   ' --maskProcessingMode ROIAUTO' ...
                   ' --ROIAutoDilateSize 3' ...
                   ' --initializeTransformMode Off' ...
                   ' --initialTransform ', ea_path_helper([volumedir, 'ct2anat.txt'])];
else
    paramset{3} = [fixparams, ...
                   ' --maskProcessingMode ROIAUTO' ...
                   ' --ROIAutoDilateSize 3'];
end
% fourth attempt..
if exist([volumedir, 'ct2anat.txt'],'file')
    paramset{4} = [fixparams, ...
                   ' --initializeTransformMode Off' ...
                   ' --initialTransform ', ea_path_helper([volumedir, 'ct2anat.txt'])];
else
    paramset{4} = fixparams;
end

basename = [fileparts(mfilename('fullpath')), filesep, 'BRAINSFit'];

if ispc
    BRAINSFit = [basename, '.exe '];
elseif isunix
    BRAINSFit = [basename, '.', computer, ' '];
end

ea_libs_helper
for trial = 1:4
    cmd = [BRAINSFit, paramset{trial}];
    if ~ispc
        status = system(['bash -c "', cmd, '"']);
    else
        status = system(cmd);
    end
    if status == 0
        fprintf(['\nBRAINSFit with parameter set ', num2str(trial), '\n']);
        break
    end
end

if ~writematout
    delete([volumedir, 'ct2anat.txt']);
end
