function output = ea_appendVolNum(niifpath, volnum)
% Append volume number to nifti file path for the use in SPM
%
% Example:
%     ea_appendVolNum('A.nii') -> 'A.nii,1'
%     ea_appendVolNum('A.nii', 1) -> 'A.nii,1'
%     ea_appendVolNum('A.nii', [1,2]) -> {'A.nii,1'; 'A.nii,2'}
%     ea_appendVolNum({'A.nii'}) -> {'A.nii,1'}
%     ea_appendVolNum({'A.nii'}, 1) -> {'A.nii,1'}
%     ea_appendVolNum({'A.nii'}, [1,2]) -> {'A.nii,1'; 'A.nii,2'}
%     ea_appendVolNum({'A.nii', 'B.nii'}, [1,2]) -> {'A.nii,1'; 'B.nii,2'}
%
% Note: Output cell will be a N*1 cell array

% Append ',1' by default
if nargin < 2
    volnum = 1;
end

% Return cell by default
outputchar = 0;
if ischar(niifpath)
    niifpath = {niifpath};
    if numel(volnum) == 1
        outputchar = 1;
    end
end

if strcmp(volnum, 'all') && numel(niifpath)==1
    volnum = 1:numel(spm_vol(niifpath{1})); % volnum = 1:end
end

% Check the sizes of 'fname' and 'volnum'
if numel(niifpath) == 1
    niifpath = repmat(niifpath, numel(volnum), 1);
elseif numel(volnum) == 1
    volnum = repmat(volnum, numel(niifpath), 1);
elseif numel(niifpath) ~= numel(volnum)
    error('''volnum'' should match the size of ''fname''!')
end

% Append volume number
output = cell(numel(niifpath), 1);
for i=1:numel(output)
    output{i} = regexprep(niifpath{i}, '\.nii$', ['.nii,', num2str(volnum(i))]);
end

% Return char as necessary
if outputchar
    output = output{1};
end
