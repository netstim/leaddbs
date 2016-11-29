function affinefile = ea_docoreg_spm(moving,fixed,cfun,doreslice,otherfiles,writeoutmat)
% Wrapper for SPM coregistration

if nargin < 3
    cfun = 'nmi';
end
if nargin < 4
    doreslice = 1;
end
if nargin < 5
    otherfiles = {''};
elseif isempty(otherfiles)  % [] or {} or ''
    otherfiles = {''};
elseif ischar(otherfiles) % single file, make it to cell string
    otherfiles = {otherfiles};
end

% Write out the transform from moving image vox to fixed image mm (also
% from fixed image vox to moving image mm)
if nargin < 6
    writeoutmat = 0;
    affinefile = {};
end

% Read the original affine matrix from the moving image (the header will be 
% overwritten after spm coregistration.
movingmat = spm_get_space(moving);

fixedmat = spm_get_space(fixed);

if doreslice
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fixed};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {moving};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = otherfiles;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = cfun;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
else
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fixed};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {moving};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = otherfiles;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = cfun;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [12 10 8 6 4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
end

spm_jobman('run',{matlabbatch});

if writeoutmat
    [~, mov] = ea_niifileparts(moving);
    [~, fix] = ea_niifileparts(fixed);
    
    % save transform: mov vox to fix mm
    spmaffine = spm_get_space(moving); % affine from mov vox to fix mm already stored in the mov header after coreg
    save([fileparts(ea_niifileparts(moving)), filesep, mov, '2', fix, '_spm.mat'], 'spmaffine', 'movingmat', 'fixedmat');
    
    % save inverse transform: fix vox to mov mm, switch fixedmat and movingmat
    spmaffine = movingmat/spm_get_space(moving)*fixedmat;
    tmp = movingmat;
    movingmat = fixedmat;
    fixedmat = tmp;
    save([fileparts(ea_niifileparts(moving)), filesep, fix, '2', mov, '_spm.mat'], 'spmaffine', 'movingmat', 'fixedmat');
    
    affinefile = {[fileparts(ea_niifileparts(moving)), filesep, mov, '2', fix, '_spm.mat'], ...
                  [fileparts(ea_niifileparts(moving)), filesep, fix, '2', mov, '_spm.mat']};
end


