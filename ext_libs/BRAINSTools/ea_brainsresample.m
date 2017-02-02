function ea_brainsresample(varargin)
% Wrapper to apply BRAINSFit transformation using BRAINSResample

fixedVolume = varargin{1};
movingVolume = varargin{2};
outputVolume = varargin{3};

volumedir = [fileparts(ea_niifileparts(movingVolume)), filesep];

if nargin >= 4
    affine = varargin{4};
else
    % determine the affine matrix to be used
    [~, mov] = ea_niifileparts(movingVolume);
    [~, fix] = ea_niifileparts(fixedVolume);
    xfm = [mov, '2', fix, '_brainsfit'];
    affine = dir([volumedir, xfm, '.h5']);

    if numel(affine) == 0
        error('Please run ea_brainsfit first before apply the transformation!');
    else
        affine = [volumedir, affine(end).name];
    end
end

% Use inverse transformation or not
if nargin >= 5
    inverse = varargin{5};
else
    inverse = 0;
end

% Available interpolation method:
% NearestNeighbor, Linear, ResampleInPlace, BSpline
% WindowedSinc, Hamming, Cosine, Welch, Lanczos, Blackman
if nargin >= 6
    interp = varargin{6};
else
    interp = 'BSpline';
end

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    BRAINSResample = [basedir, 'BRAINSResample.exe'];
else
    BRAINSResample = [basedir, 'BRAINSResample.', computer('arch')];
end

cmd = [BRAINSResample, ...
       ' --referenceVolume ', ea_path_helper(fixedVolume), ...
       ' --inputVolume ', ea_path_helper(movingVolume), ...
       ' --outputVolume ', ea_path_helper(outputVolume), ...
       ' --warpTransform ', ea_path_helper(affine), ...
       ' --interpolationMode ', interp, ...
       ' --pixelType float'];

if inverse
    cmd = [cmd, ' --inverseTransform'];
end

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
