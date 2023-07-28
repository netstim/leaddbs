function header = ea_fslhd(input, xmlarg)
% Wrapper for fslhd

% Use XML-style NIFTI header or not
if nargin < 2
    xmlarg = '';
elseif strcmp(xmlarg, 'x')
    xmlarg = '-x';
end

if ~isfile(input)
    error('%s not found!', input);
end

basedir = [fileparts(mfilename('fullpath')), filesep];
FSLHD = ea_getExec([basedir, 'fslhd'], escapePath = 1);


cmd = [FSLHD, ' ', xmlarg, ' ', ea_path_helper(input), ];
[status, cmdout] = ea_runcmd(cmd);

if status ~= 0
    error('%s', strip(cmdout));
end

% Trim string
if strcmp(xmlarg, '-x')
    cmdout = regexprep(cmdout, '(<nifti_image\n|  |\n/>\n)', '');
    cmdout = regexprep(cmdout, '= ''', '= ');
    cmdout = regexprep(cmdout(1:end-1), '''\n', '\n');
else
    cmdout  = strrep(cmdout(1:end-1), 'size of header', 'sizeof_header');
    cmdout  = strrep(cmdout, 'xyz:', 'xyz');
    cmdout = regexprep(cmdout, '[\t]+', ' = ');
end

cmdout = cellfun(@(x) strsplit(x, ' = '), strsplit(cmdout, ea_newline), 'Uni', 0)';

% Construct header
header = struct;
for i=1:length(cmdout)
    if length(cmdout{i})==1
        continue % shell init issue on macOS, skip this line
    else
        if isempty(cmdout{i}{2})
            eval(['header.', cmdout{i}{1}, ' = '''';']);
        elseif ~isempty(regexp(cmdout{i}{2}, '[0-9]', 'once')) && ...
                isempty(regexp(cmdout{i}{2}, '[^0-9.+-e ]', 'once'))
            try
                eval(['header.', cmdout{i}{1}, ' = [', cmdout{i}{2}, '];']);
            catch
                eval(['header.', cmdout{i}{1}, ' = ''', cmdout{i}{2}, ''';']);
            end
        else
            eval(['header.', cmdout{i}{1}, ' = ''', cmdout{i}{2}, ''';']);
        end
    end
end

% Reshape sto_xyz_matrix to 4x4 size
if isfield(header, 'sto_xyz_matrix')
    header.sto_xyz_matrix = reshape(header.sto_xyz_matrix,4,4)';
end
