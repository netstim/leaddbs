function pointsw = ea_antspy_apply_transforms_to_points(points, transform, useinverse)
% Wrapper for ANTsPy's apply_transforms_to_points function

arguments
    points       {mustBeNumeric}                   % N*3 Input coordinates
    transform    {mustBeText}                      % Transform or Lead-DBS patient folder
    useinverse   {mustBeNumericOrLogical} = false  % Use inverse transform or not
end

if size(points,2) ~= 3
    ea_error('Input must have the size of N*3!', showdlg=0, simpleStack=1);
end

% Lead-DBS subj folder instead of transfom provided
if isfolder(transform)
    options = ea_getptopts(transform);

    json = loadjson(options.subj.norm.log.method);
    if contains(json.method, 'affine', 'IgnoreCase', true)
        % Three-step affine normalization (Schonecker 2009) used
        warpSuffix = 'ants.mat';
    else
        warpSuffix = 'ants.nii.gz';
    end

    if useinverse
        transform = [options.subj.norm.transform.inverseBaseName, warpSuffix];
    else
        transform = [options.subj.norm.transform.forwardBaseName, warpSuffix];
    end

    useinverse = false;
end

if ~isfile(transform)
    ea_error('Transformation not found!', showdlg=0, simpleStack=1);
end

% Check Lead-DBS conda env for ANTsPy
condaenv = ea_conda_env('Lead-DBS');
if ~condaenv.is_created
    ea_cprintf('*Comments', 'Initializing Lead-DBS conda environment...\n')
    condaenv.create;
    ea_cprintf('*Comments', 'Lead-DBS conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('*Comments', 'Updating Lead-DBS conda environment...\n')
    condaenv.update;
    ea_cprintf('*Comments', 'Lead-DBS conda environment initialized.\n')
end

% Set pyenv
pe = pyenv;
restoreENV = 0;
if pe.Executable ~= condaenv.python
    if pe.Version ~= ""
        restoreENV = 1;
    end
    pyenv('Version', condaenv.python);
end

% Prepare input data frame for ANTsPy
df = py.pandas.DataFrame(py.numpy.array(points).reshape(py.int(-1),py.int(3)), pyargs('columns', {'x', 'y', 'z'}));

% Apply transform to points
pointsw = py.ants.apply_transforms_to_points(py.int(3), df, transform, py.list({useinverse}));
pointsw = cast(double(pointsw.values), class(points));

% Restore pyenv
if restoreENV
    pyenv('Version', pe.Executable);
end
