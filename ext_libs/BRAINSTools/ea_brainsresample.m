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

% Available interpolation method:
% NearestNeighbor, Linear, ResampleInPlace, BSpline
% WindowedSinc, Hamming, Cosine, Welch, Lanczos, Blackman
if nargin >= 5
    interp = varargin{5};
else
    interp = 'BSpline';
end

% Use inverse transformation or not
if nargin >= 6
    inverse = varargin{6};
else
    inverse = 0;
end

basedir = [fileparts(mfilename('fullpath')), filesep];
BRAINSResample = ea_getExec([basedir, 'BRAINSResample'], escapePath = 1);


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

ea_runcmd(cmd);
