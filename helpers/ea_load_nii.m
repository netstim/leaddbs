function nii=ea_load_nii(fname, mode)
% simple nifti reader using SPM.

% try to combine ea_load_nii and load_nii_proxi here, temporarily solution
% will have a unified interface for nii/nii.gz reading in the future
% by default, use the original ea_load_nii way
% if mode=='simple', use the load_nii_proxi way
if nargin < 2
    mode = 0; 
end

% need to consider the spm_vol('image.nii,1') case
fname=[ea_niigz(fname), fname(strfind(fname, ','):end)];

if regexp(fname, '\.nii.gz$', 'once') % 'image.nii.gz'
    wasgzip=1;
    gunzip(fname);
    fname=fname(1:end-3);
elseif regexp(fname, '\.nii.gz,\d+$', 'once') % 'image.nii.gz,1'
    wasgzip=1;
    [pth, ~, ~, slice] = ea_niifileparts(fname);
    gunzip([pth, '.nii.gz']);
    fname = [pth, '.nii', slice];
else
    wasgzip=0;
end

if mode
    V=spm_vol(fname);
    X=spm_read_vols(V);
    nii.img=X;
    nii.hdr=V;
    nii.hdr.dime.pixdim=ea_detvoxsize(V(1).mat);
else
    nii=spm_vol(fname);
    img=spm_read_vols(nii);
    if length(nii)>1 % multi volume;
        anii=nii(1);
        anii.img=img;
        nii=anii;
        clear anii
    else
        nii.img=img;
    end
    try
        nii.hdr.dime.pixdim=ea_detvoxsize(nii(1).mat);
    catch
        fprintf('Can not determine the voxel size..\n')
    end
end

if wasgzip
    delete([ea_niifileparts(fname), '.nii.gz']); % since gunzip makes a copy of the zipped file.
end
