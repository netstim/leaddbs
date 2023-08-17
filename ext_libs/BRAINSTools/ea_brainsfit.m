function affinefile = ea_brainsfit(varargin)
% Wrapper for BRAINSFit

fixedVolume = varargin{1};
movingVolume = varargin{2};
outputVolume = varargin{3};

if nargin >= 4
    writeoutmat = varargin{4};
else
    writeoutmat = 1;
end

if nargin >= 5
    if isempty(varargin{5}) % [] or {} or ''
        otherfiles = {};
    elseif ischar(varargin{5}) % single file, make it to cell string
        otherfiles = varargin(5);
    else % cell string
        otherfiles = varargin{5};
    end
else
    otherfiles = {};
end

volumedir = [fileparts(ea_niifileparts(outputVolume)), filesep];

fixedVolume = ea_niigz(fixedVolume);
movingVolume = ea_niigz(movingVolume);
outputVolume = ea_niigz(outputVolume);

% name of the output transformation
[~, mov] = ea_niifileparts(movingVolume);
[~, fix] = ea_niifileparts(fixedVolume);
xfm = [mov, '2', fix, '_brainsfit'];

fixparams = [' --fixedVolume ' , ea_path_helper(fixedVolume), ...
             ' --movingVolume ', ea_path_helper(movingVolume), ...
             ' --outputVolume ', ea_path_helper(outputVolume), ...
             ' --useRigid --useAffine' ...
             ' --samplingPercentage 0.005' ...
             ' --removeIntensityOutliers 0.005' ...
             ' --interpolationMode Linear' ...
             ' --outputTransform ', ea_path_helper([volumedir, xfm, '.h5'])];

% first attempt...
paramset{1} = [fixparams, ...
               ' --initializeTransformMode useGeometryAlign' ...
               ' --maskProcessingMode ROIAUTO' ...
               ' --ROIAutoDilateSize 3'];

% second attempt...
paramset{2} = [fixparams, ...
               ' --initializeTransformMode useGeometryAlign'];

% third attempt...
if exist([volumedir, xfm, '.h5'],'file')
    paramset{3} = [fixparams, ...
                   ' --maskProcessingMode ROIAUTO' ...
                   ' --ROIAutoDilateSize 3' ...
                   ' --initializeTransformMode Off' ...
                   ' --initialTransform ', ea_path_helper([volumedir, xfm, '.h5'])];
else
    paramset{3} = [fixparams, ...
                   ' --maskProcessingMode ROIAUTO' ...
                   ' --ROIAutoDilateSize 3'];
end

% fourth attempt..
if exist([volumedir, xfm, '.h5'],'file')
    paramset{4} = [fixparams, ...
                   ' --initializeTransformMode Off' ...
                   ' --initialTransform ', ea_path_helper([volumedir, xfm, '.h5'])];
else
    paramset{4} = fixparams;
end

basename = [fileparts(mfilename('fullpath')), filesep, 'BRAINSFit'];
BRAINSFit = ea_getExec(basename, escapePath = 1);


ea_libs_helper
for trial = 1:4
    cmd = [BRAINSFit, ' ', paramset{trial}];
    status = ea_runcmd(cmd);

    if status == 0
        fprintf(['\nBRAINSFit with parameter set ', num2str(trial), '\n']);
        break
    end
end

if ~isempty(otherfiles)
    for fi = 1:numel(otherfiles)
        ea_brainsresample(fixedVolume, otherfiles{fi}, otherfiles{fi}, [volumedir, xfm, '.h5']);
    end
end

if ~writeoutmat
    ea_delete([volumedir, xfm, '.h5']);
    ea_delete([volumedir, xfm, '_Inverse.h5']);
    affinefile = {};
else
    % Convert h5 to mat
    AffineTransform_double_3_3 = h5read([volumedir, xfm, '.h5'], '/TransformGroup/0/TransformParameters');
    fixed = h5read([volumedir, xfm, '.h5'], '/TransformGroup/0/TransformFixedParameters');
    save([volumedir, xfm, '.mat'], 'AffineTransform_double_3_3', 'fixed');
    delete([volumedir, xfm, '.h5']);

    % Convert h5 to mat
    AffineTransform_double_3_3 = h5read([volumedir, xfm, '_Inverse.h5'], '/TransformGroup/0/TransformParameters');
    fixed = h5read([volumedir, xfm, '_Inverse.h5'], '/TransformGroup/0/TransformFixedParameters');
    save([volumedir, xfm, '_Inverse.mat'], 'AffineTransform_double_3_3', 'fixed');
    delete([volumedir, xfm, '_Inverse.h5']);

    invxfm = [fix, '2', mov, '_brainsfit'];
    movefile([volumedir, xfm, '_Inverse.mat'], [volumedir, invxfm, '.mat']);

    affinefile = {[volumedir, xfm, '.mat']
                  [volumedir, invxfm, '.mat']};
end

fprintf('\nBRAINSFit LINEAR registration done.\n');

%% add methods dump:
cits={
    'Johnson, H., Harris, G., & Williams, K. (2007). BRAINSFit: mutual information rigid registrations of whole-brain 3D images, using the insight toolkit. Insight J.'
};

ea_methods(volumedir,[mov,' was linearly co-registered to ',fix,' using BRAINSFit software (Johnson 2007; https://www.nitrc.org/projects/multimodereg/)'],...
    cits);
