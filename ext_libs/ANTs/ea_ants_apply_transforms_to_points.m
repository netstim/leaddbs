function output = ea_ants_apply_transforms_to_points(varargin)
% Wrapper for antsApplyTransformsToPoints used for both linear
% transformation and non-linear transformation

subjDir = varargin{1};
input = varargin{2};
useinverse = varargin{3};

options = ea_getptopts(subjDir);

if nargin>3 % Transformation file explicitly specified
    transform = varargin{4};
    tstring = [' --transform [', ea_path_helper(transform), ',',num2str(useinverse),']']; % [transformFileName,useInverse]
else
    if useinverse
        transform = [options.subj.norm.transform.inverseBaseName, 'ants.nii.gz'];
    else
        transform = [options.subj.norm.transform.forwardBaseName, 'ants.nii.gz'];
    end

    if isfile(transform)
        tstring = [' -t [', ea_path_helper(transform), ',0]'];
    else
        error('Transformation file not found! Please rerun normalization.');
    end
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransformsToPoints = ea_path_helper([basedir, 'antsApplyTransformsToPoints.exe']);
else
    applyTransformsToPoints = ea_path_helper([basedir, 'antsApplyTransformsToPoints.', computer('arch')]);
end

tmpdir = fullfile(subjDir, 'tmp');
ea_mkdir(tmpdir);
uuid = ea_generate_uuid;

cmd = [applyTransformsToPoints, ...
    ' --dimensionality 3' ...   % dimensionality
    ' --precision 0' ...    % single precision
    ' --input ', ea_path_helper([tmpdir, filesep, 'tmpin_', uuid, '.csv']) ...  % input csv file with x,y,z,t (at least) as the column header
    ' --output ', ea_path_helper([tmpdir, filesep, 'tmpout_', uuid, '.csv']) ...    % warped output csv file
    tstring];

ea_writecsv([tmpdir, filesep, 'tmpin_', uuid, '.csv'], input);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

output = ea_readcsv([directory,'tmpout_',uuid,'.csv']);
rmdir(tmpdir, 's');


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
