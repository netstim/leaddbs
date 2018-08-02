function ea_spm_coreg_applytransform(varargin)
% Wrapper to apply SPM coregistration transformation

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

directory = fileparts(ea_niifileparts(movingimage));

if nargin >= 4
    affine = varargin{4};
else
    % determine the affine matrix to be used
    [~, mov] = ea_niifileparts(movingimage);
    [~, fix] = ea_niifileparts(fixedimage);
    affine = [directory, filesep, mov, '2', fix, '_spm.mat'];
    if ~exist(affine, 'file')
        error('No proper SPM transformation found!');
    end
end

if nargin >= 5
    interp = varargin{5};
else
    interp = 1; % Linear interpolition by default
end

% set header affine matrix
load(affine, 'spmaffine');
nii = ea_load_nii(movingimage);
nii.mat = spmaffine;
nii.fname = outputimage;
ea_write_nii(nii);

% reslice
matlabbatch{1}.spm.spatial.coreg.write.ref = {fixedimage};
matlabbatch{1}.spm.spatial.coreg.write.source = {outputimage};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = '';
spm_jobman('run',{matlabbatch});
clear matlabbatch
 