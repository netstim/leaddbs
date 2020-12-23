function outcoords = ea_fsl_img2imgcoord(varargin)
% Wrapper for FSL img2imgcoord
% Input coords should have N*3 size

incoords = varargin{1};
src = varargin{2};
dest = varargin{3};
transform = varargin{4};

% Apply warp field or affine transform
if strcmp(varargin{5}, 'l')
    type = ' -xfm';
elseif strcmp(varargin{5}, 'n')
    type = ' -inv -warp'; % inverse warp field supplied
end

% Input coords in mm by default
if nargin >= 6
    space = varargin{6};
else
    space = 'mm';
end

if nargin >= 9
    premat = varargin{9};
else
    premat = '';
end

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    IMG2IMGCOORD = ea_path_helper([basedir, 'img2imgcoord.exe']);
else
    IMG2IMGCOORD = [basedir, 'img2imgcoord.', computer('arch')];
end

cmd = [IMG2IMGCOORD, ...
       ' -src ' ea_path_helper(src), ...
       ' -dest ' ea_path_helper(dest), ...
       type, ' ', ea_path_helper(transform)];

if strcmp(space, 'mm')
    cmd = [cmd, ' -mm'];
else
    cmd = [cmd, ' -vox'];
end

if ~isempty(premat)
    cmd = [cmd, ' -premat ', ea_path_helper(premat)];
end

% Input coords is var rather than a file
if ~ischar(incoords)
    % convert column vector to make it of the size N*3
    if size(incoords, 2) ~= 3
        incoords = incoords';
    end
    % Write temporary coords file
    uuid=ea_generate_uuid;
    directory = [fileparts(ea_niifileparts(src)), filesep];
    fid=fopen([directory,'tmpin_',uuid,'.csv'],'w');
    fprintf(fid,'%.9f %.9f %.9f\n', incoords'); % transpose needed for 'fprintf': matrix column to file row.
    fclose(fid);
    incoords = [directory,'tmpin_',uuid,'.csv'];
end

cmd = [cmd, ' ', ea_path_helper(incoords)];

setenv('FSLOUTPUTTYPE', 'NIFTI');
if ~ispc
    [status, cmdout] = system(['bash -c "', cmd, '"']);
else
    [status, cmdout] = system(cmd);
end

if status == 0
    outcoords = cell2mat(textscan(cmdout, '%f %f %f', 'HeaderLines', 1));
else
    error(['Coords mapping failed:\n', cmdout]);
end

% Delete tmp file
if exist('uuid', 'var')
    delete(incoords);
end
