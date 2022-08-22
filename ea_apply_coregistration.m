function ea_apply_coregistration(varargin)
% Wrapper to apply linear transformation
%
% The 4th parameter (transformation) can be:
%     1. Not specified, the transformation will be automatically detected
%        if there is only one transformation file in the folder.
%     2. String ('SPM', 'ANTs', 'FSL' and 'BRAINSFIT'), the transformation
%        will be chosen according to the type specified in the string.
%     3. String (path to the transformation file).

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

[~, mov] = ea_niifileparts(movingimage);
[~, fix] = ea_niifileparts(fixedimage);

volumedir = [fileparts(ea_niifileparts(movingimage)), filesep];

if nargin < 4
    % determine the transformation to be used
    xfm = [mov, '2', fix, '_\w+\.(mat|h5)$'];
    transform = ea_regexpdir(volumedir, xfm, 0);

    if numel(transform) == 0
        error(['Transformation not found! Please run coregistration ' ...
               'first before applying the transformation!']);
    elseif numel(transform) > 1
        error(['Multiple transformations found! Please explicitly ' ...
               'specify the one to be used in the 4th parameter.']);
    else % Only one transformation found
        transform = transform{1};
    end
elseif regexp(varargin{4}, '^[a-zA-Z]+( +[a-zA-Z]+)*$')
% elseif isfile(varargin{4})
    transformType = lower(varargin{4});
    if regexp(transformType, '^fsl') % Determine FSL transformation name
        transformType = regexp(transformType, '(?<=^fsl )(.+)$', 'match', 'once');
    end

    xfm = [mov, '2', fix, '_', transformType, '\d*\.(mat|h5)$'];
    transform = ea_regexpdir(volumedir, xfm, 0);

    if numel(transform) == 0
        error(['Specified transformation not found! Please run ' ...
               'coregistration first before applying the transformation!']);
    else
        if numel(transform) > 1
            warning(['Multiple transformations of the same type found! ' ...
                     'Will use the last one:\n%s'], transform{end});
        end
        transform = transform{end};
    end
else
    transform = varargin{4};
    if ~exist(transform, 'file')
        error(['Specified transformation not found! Please run ' ...
               'coregistration first before applying the transformation!']);
    end
end

transformType = upper(regexp(transform, '(?<=desc-)\w+(?=\.mat)', 'match', 'once'));

% nn: NearestNeighbor
% lbl: Label (for multi-label image, ANTs uses 'GenericLabel' rather than 'NearestNeighbor')
% ln: Linear
% spl: Spline
if nargin < 5
    interp = 'ln'; % Linear interpolation by default
else
    interp = varargin{5};
end

% Determine interpolation type
switch transformType
    case 'SPM'
        % 0: Nearest neighbour
        % 1: Trilinear
        % 2: 2nd Degree B-Spline
        % 3: 3nd Degree B-Spline
        % 4: 4nd Degree B-Spline
        if ischar(interp)
            switch lower(interp)
                case {'nn', 'nearestneighbor', 'lbl', 'label'}
                    interp = 0;
                case {'ln', 'linear'}
                    interp = 1;
                case {'spl', 'spline'}
                    interp = 4; % default degree in SPM coregistration
            end
        end
    case 'ANTS'
        % Linear, NearestNeighbor, MultiLabel, Gaussian, BSpline
        % CosineWindowedSinc, WelchWindowedSinc, HammingWindowedSinc, LanczosWindowedSinc
        % GenericLabel (Recommanded for label image)
        switch lower(interp)
            case {'nn', 'nearestneighbor'}
                interp = 'NearestNeighbor';
            case {'lbl', 'label'}
                interp = 'GenericLabel';
            case {'ln', 'linear'}
                interp = 'Linear';
            case {'spl', 'spline'}
                interp = 'BSpline';
        end
    case {'FLIRT','FLIRTBBR', 'BBR', 'FSL'}
        % trilinear, nearestneighbour, sinc, spline
        switch lower(interp)
            case {'nn', 'nearestneighbor', 'lbl', 'label'}
                interp = 'nearestneighbour';
            case {'ln', 'linear'}
                interp = 'trilinear';
            case {'spl', 'spline'}
                interp = 'spline';
            otherwise
                error('Unrecognized interpolation type: ''%s''!', interp);
        end
    case 'BRAINSFIT'
        % NearestNeighbor, Linear, ResampleInPlace, BSpline
        % WindowedSinc, Hamming, Cosine, Welch, Lanczos, Blackman
        switch lower(interp)
            case {'nn', 'nearestneighbor', 'lbl', 'label'}
                interp = 'NearestNeighbor';
            case {'ln', 'linear'}
                interp = 'Linear';
            case {'spl', 'spline'}
                interp = 'BSpline';
            otherwise
                error('Unrecognized interpolation type: ''%s''!', interp);
        end
    otherwise
        error('Unrecognized transformation type: ''%s''!', transformType);
end

% Apply transform
switch transformType
    case 'SPM'
        ea_spm_apply_coregistration(fixedimage, movingimage, outputimage, ...
            transform, interp)
    case 'ANTS'
        ea_ants_apply_transforms([], movingimage, outputimage, ...
            0, fixedimage, transform, interp)
    case 'FSL'
        ea_fsl_apply_coregistration(fixedimage, movingimage, outputimage, ...
            transform, interp)
    case 'FLIRT'
        ea_fsl_apply_coregistration(fixedimage, movingimage, outputimage, ...
            transform, interp)
    case 'BRAINSFIT'
        ea_brainsresample(fixedimage, movingimage, outputimage, ...
            transform, interp)
    otherwise
        error('Unrecognized transformation type: ''%s''!', transformType);
end
