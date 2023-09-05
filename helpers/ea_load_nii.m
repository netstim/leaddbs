function nii = ea_load_nii(fname)
% Simple nifti reader using SPM.
%
% For multi-vol image, the data of all the vols will be read in, but only
% the header of the first vol will be kept in the returned structure.
%
% Support fname:
%     'PATH/TO/image.nii'
%     'PATH/TO/image.nii,1' (volume 1)
%     'PATH/TO/image.nii.gz'
%     'PATH/TO/image.nii.gz,1'
%     'PATH/TO/image'
%     'PATH/TO/image,1'

% input file is MAT file
if endsWith(fname, '.mat')
    nii = load(fname);
    return
end

% need to consider the case like: spm_vol('image.nii,1')
fname = [ea_niigz(fname), fname(strfind(fname, ','):end)];

 % 'image.nii.gz' or 'image.nii.gz,1'
if regexp(fname, '\.nii.gz(,\d+)?$', 'once')
    wasgzip = 1;
    gunzip(regexprep(fname, ',\d+$', ''));
    fname = regexprep(fname, '\.gz(,\d+)?$', '$1');
else
    wasgzip = 0;
end

% Read header and img

nii = spm_vol(fname);
img = spm_read_vols(nii);
volnum = numel(nii);
nii = nii(1); % only keep the header of the first vol for multi-vol image
nii.img = img; % set image data
nii.voxsize = ea_detvoxsize(nii.mat); % set voxsize
nii.volnum = volnum; % set number of volumes

if wasgzip
    delete(regexprep(fname, ',\d+$', '')); % since gunzip makes a copy of the zipped file.
end
