function header = ea_fslhd(input, xmlarg)
% Wrapper for fslhd

% Use XML-style NIFTI header or not
if nargin < 2
    xmlarg = '';
elseif strcmp(xmlarg, 'x')
    xmlarg = '-x';
end

input = ea_path_helper(input);

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    FSLHD = ea_path_helper([basedir, 'fslhd.exe']);
else
    FSLHD = [basedir, 'fslhd.', computer('arch')];
end

cmd = [FSLHD, ' ', xmlarg, ' ', input, ];

if ~ispc
    [~, cmdout] = system(['bash -c "', cmd, '"']);
else
    [~, cmdout] = system(cmd);
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

% Reshape sto_xyz_matrix to 4x4 size
if isfield(header, 'sto_xyz_matrix')
    header.sto_xyz_matrix = reshape(header.sto_xyz_matrix,4,4)';
end
