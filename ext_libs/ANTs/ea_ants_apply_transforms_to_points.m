function output = ea_ants_apply_transforms_to_points(varargin)
% Wrapper for antsApplyTransformsToPoints used for both linear
% transformation and non-linear transformation

subjDir = varargin{1};
input = varargin{2};
useinverse = varargin{3};

if nargin>3 % Transformation file explicitly specified
    transform = varargin{4};
    tstring = [' --transform [', ea_path_helper(transform), ',',num2str(useinverse),']']; % [transformFileName,useInverse]
else
    options = ea_getptopts(subjDir);

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

    if isfile(transform)
        tstring = [' -t [', ea_path_helper(transform), ',0]'];
    else
        error('Transformation file not found! Please rerun normalization.');
    end
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

applyTransformsToPoints = ea_getExec([basedir, 'antsApplyTransformsToPoints'], escapePath = 1);


uuid = ea_generate_uuid;

input_file = [tempdir, 'tmpin_', uuid, '.csv'];
output_file = [tempdir, 'tmpout_', uuid, '.csv'];

cmd = [applyTransformsToPoints, ...
    ' --dimensionality 3' ...   % dimensionality
    ' --precision 0' ...    % single precision
    ' --input ',  ea_path_helper(input_file) ...  % input csv file with x,y,z,t (at least) as the column header
    ' --output ', ea_path_helper(output_file) ...    % warped output csv file
    tstring];

ea_writecsv(input_file, input);

ea_runcmd(cmd);

output = ea_readcsv(output_file);
ea_delete(input_file);
ea_delete(output_file)


function coord = ea_readcsv(pth)
fid = fopen(pth);
C = textscan(fid,'%f %f %f %f','commentStyle', '#','delimiter', ',','Headerlines',1);
fclose(fid);
coord = cell2mat(C(1:3));


function ea_writecsv(pth,input)
fid = fopen(pth,'w');
try
    fprintf(fid,'x,y,z,t \n');
catch
    ea_error(['Cannot open file for writing at ',pth,'.']);
end
fprintf(fid,'%.9f, %.9f, %.9f, 0\n',input'); % transpose needed for 'fprintf': matrix column to file row
fclose(fid);
