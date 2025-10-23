function output = ea_prependFilename(filepath, prefix, suffix)
% EA_PREPENDFILENAME  Add prefix and/or suffix to filename (not path)
%
% For BIDS paths like 'preprocessing/func/sub-XXX_bold.nii':
%   ea_prependFilename('preprocessing/func/sub-XXX_bold.nii', 'r')
%     -> 'preprocessing/func/rsub-XXX_bold.nii'
%
% For classic paths like 'rest.nii':
%   ea_prependFilename('rest.nii', 'r')
%     -> 'rrest.nii'
%
% With suffix:
%   ea_prependFilename('preprocessing/func/sub-XXX_bold.nii', '', '_anat_t2.nii')
%     -> 'preprocessing/func/sub-XXX_bold_anat_t2.nii'

if nargin < 2
    prefix = '';
end

if nargin < 3
    suffix = '';
end

% Split path and filename
[fpath, fname, ext] = fileparts(filepath);

% Add prefix to filename
output = fullfile(fpath, [prefix, fname, suffix, ext]);

% Remove leading filesep if original didn't have path
if isempty(fpath) && startsWith(output, filesep)
    output = output(2:end);
end
end

