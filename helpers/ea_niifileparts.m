function [trimpath, basename, ext, vol] = ea_niifileparts(niifile)
% Return nifti file path and name without .nii/.nii.gz extension.
% Useful for FSL tools.
% For example, input '/PATH/TO/image.nii.gz,1' will return 
% ['/PATH/TO/image', 'image', '.nii.gz', ',1']

% TP - fixing bug:
% on Windows PC, paths are encapsulated in double quotes (see
% ea_path_helper.m). The regexp fails since end of line does not occur
% after extension. My solution is not very elegant - perhaps there is a
% better fix.

is_quote = 0;
if ispc
    if niifile(1) == '"'
        niifile = niifile(2:end-1); % remove double quotes and add after regex
        is_quote = 1;
    end
end

if regexp(niifile, '\.nii$', 'once') % '/PATH/TO/image.nii'
    trimpath = niifile(1:end-4);
    ext = '.nii';
    vol = '';
elseif ~isempty(regexp(niifile, '(?<=\.nii),\d+$', 'match')) % '/PATH/TO/image.nii,1'
    trimpath = niifile(1:regexp(niifile, '\.nii,\d+$', 'once')-1);
    ext = '.nii';
    vol = regexp(niifile, '(?<=\.nii),\d+$', 'match');
    vol = vol{:};
elseif regexp(niifile, '\.nii.gz$', 'once') % '/PATH/TO/image.nii.gz'
    trimpath = niifile(1:end-7);
    ext = '.nii.gz';
    vol = '';
elseif ~isempty(regexp(niifile, '(?<=\.nii.gz),\d+$', 'match')) % '/PATH/TO/image.nii.gz,1'
    trimpath = niifile(1:regexp(niifile, '\.nii.gz,\d+$', 'once')-1);
    ext = '.nii.gz';
    vol = regexp(niifile, '(?<=\.nii.gz),\d+$', 'match');
    vol = vol{:};
elseif ~isempty(regexp(niifile, ',\d+$', 'match')) % '/PATH/TO/image,1'
    trimpath = niifile(1:regexp(niifile, ',\d+$', 'once')-1);
    ext = '';
    vol = regexp(niifile, ',\d+$', 'match');
    vol = vol{:};
else % '/PATH/TO/image'
    trimpath = niifile;
    ext = '';
    vol = '';
end

if isempty(fileparts(niifile))
    trimpath = ['.', filesep, trimpath];
end

if is_quote
    trimpath = ['"',trimpath,'"']; % add double quotes if originally present
end

basename = regexprep(trimpath, ['.*\', filesep], '');
