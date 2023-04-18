function fn = ea_rmext(fn)
% Remove .nii/.nii.gz extension form file name[s]

% if iscell(fn)
%     [~, fn] = cellfun(@fileparts, fn, 'UniformOutput', 0); % loose .gz
%     [~, fn] = cellfun(@fileparts, fn, 'UniformOutput', 0); % loose .nii
% else
%     [~, fn] = fileparts(fn); % loose .gz
%     [~, fn] = fileparts(fn); % loose .nii
% end

if iscell(fn)
    fn = cellfun(@(x) regexp(x, '(.*)(?=\.nii(\.gz)?$)', 'match', 'once'), fn, 'UniformOutput', 0);
else
    fn = regexp(fn, '(.*)(?=\.nii(\.gz)?$)', 'match', 'once');
end

