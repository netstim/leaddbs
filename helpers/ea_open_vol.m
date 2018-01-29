function V = ea_open_vol(fname)
% Enhanced wrapper for 'spm_vol'.
%
% For multi-vol image, only the header of the first vol will be kept in the
% returned structure.
%
% Support fname:
%     'PATH/TO/image.nii'
%     'PATH/TO/image.nii,1' (volume 1)
%     'PATH/TO/image.nii.gz'
%     'PATH/TO/image.nii.gz,1'
%     'PATH/TO/image'
%     'PATH/TO/image,1'

% need to consider the case like: spm_vol('image.nii,1')
fname = [ea_niigz(fname), fname(strfind(fname, ','):end)];

if regexp(fname, '\.nii.gz$', 'once') % 'image.nii.gz'
    wasgzip = 1;
    gunzip(fname);
    fname = fname(1:end-3);
elseif regexp(fname, '\.nii.gz,\d+$', 'once') % 'image.nii.gz,1'
    wasgzip = 1;
    [pth, ~, ~, vol]  =  ea_niifileparts(fname);
    gunzip([pth, '.nii.gz']);
    fname  =  [pth, '.nii', vol];
else
    wasgzip = 0;
end

% Read header
V = spm_vol(fname);
volnum = numel(V);
V = V(1); % only keep the header of the first vol for multi-vol image
V.voxsize = ea_detvoxsize(V.mat); % set voxsize
V.volnum = volnum; % set number of volumes

if wasgzip
    delete(fname); % since gunzip makes a copy of the zipped file.
end
