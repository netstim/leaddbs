function affinefile = ea_spm_coreg(options,moving,fixed,costfun,doreslice,otherfiles,writeoutmat,interp)
% Wrapper for SPM coregistration
%
% SPM modifies the header of the moving image when doing coregistration
% (for both "Estimate" and "Estimate & Reslice"). This wrapper here will
% only modify the header if 'doreslice' is FALSE (in-place "Estimate"). If
% 'doreslice' is TRUE ("Estimate & Reslice"), the moving image will be left
% untouched and a new coregistered image will be generated with 'r' prefix.

if ~exist('costfun','var')
    costfun = 'nmi';
elseif isempty(costfun)
    costfun = 'nmi';
end

if ~exist('doreslice','var')
    doreslice = 1;
elseif isempty(doreslice)
    doreslice = 1;
end

if ~exist('otherfiles','var')
    otherfiles = {''};
elseif isempty(otherfiles)  % [] or {} or ''
    otherfiles = {''};
elseif ischar(otherfiles) % single file, make it to cell string
    otherfiles = {otherfiles};
end

if ~exist('interp','var')
    interp = 4;
elseif isempty(interp)
    interp = 4;
end

% Write out the transform from moving image vox to fixed image mm (also
% from fixed image vox to moving image mm)
if ~exist('writeoutmat','var')
    writeoutmat = 0;
elseif isempty(writeoutmat)
    writeoutmat = 0;
end

% Read the original affine matrix from the moving image (the header will be
% overwritten after spm coregistration.
movingmat = spm_get_space(moving);
fixedmat = spm_get_space(fixed);

if doreslice
    % backup moving image
    movingbackup = ea_niifileparts(moving);
    copyfile(regexprep(moving, ',\d+$', ''), movingbackup);

    % backup other files
    otherfilesbackup = cell(size(otherfiles));
    for i=1:length(otherfiles)
        if ~isempty(otherfiles{i})
            otherfilesbackup{i} = ea_niifileparts(otherfiles{i});
            copyfile(regexprep(otherfiles{i}, ',\d+$', ''), otherfilesbackup{i});
        end
    end

    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ea_appendVolNum(fixed)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {ea_appendVolNum(moving)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = ea_appendVolNum(otherfiles);
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = costfun;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [8 4 2]; %[4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = interp;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',{matlabbatch});
else
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fixed};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {moving};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = otherfiles;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfun;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [12 10 8 6 4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run',{matlabbatch});
end

ea_delete([fileparts(ea_niifileparts(moving)), filesep, 'spm*.ps']);

[~, mov] = ea_niifileparts(moving);
[~, fix] = ea_niifileparts(fixed);

if ~writeoutmat
    affinefile = {};
else
    % save transform: mov vox to fix mm
    spmaffine = spm_get_space(moving); % affine from mov vox to fix mm already stored in the mov header after coreg
    tmat = inv(movingmat/spmaffine); % matrix from moving mm to fix mm
    save([fileparts(ea_niifileparts(moving)), filesep, mov, '2', fix, '_spm.mat'], 'spmaffine', 'movingmat', 'fixedmat', 'tmat');

    % save inverse transform: fix vox to mov mm, switch fixedmat and movingmat
    spmaffine = movingmat/spm_get_space(moving)*fixedmat;
    tmp = movingmat;
    movingmat = fixedmat;
    fixedmat = tmp;
    tmat = inv(movingmat/spmaffine); % matrix from fixed mm to moving mm
    save([fileparts(ea_niifileparts(moving)), filesep, fix, '2', mov, '_spm.mat'], 'spmaffine', 'movingmat', 'fixedmat', 'tmat');

    affinefile = {[fileparts(ea_niifileparts(moving)), filesep, mov, '2', fix, '_spm.mat']
        [fileparts(ea_niifileparts(moving)), filesep, fix, '2', mov, '_spm.mat']};
end

if doreslice
    % restore moving image
    movefile(movingbackup, regexprep(moving, ',\d+$', ''));

    % restore other files
    for i=1:length(otherfilesbackup)
        if ~isempty(otherfilesbackup{i})
            movefile(otherfilesbackup{i}, regexprep(otherfiles{i}, ',\d+$', ''));
        end
    end
end

%% add methods dump:
cits={
    'Friston, K. J., Ashburner, J. T., Kiebel, S. J., Nichols, T. E., & Penny, W. D. (2007). Statistical Parametric Mapping: The Analysis of Functional Brain Images. Academic Press.'
};

ea_methods(options,[mov,' was linearly co-registered to ',fix,' using ',spm('ver'),' (Friston 2007; http://www.fil.ion.ucl.ac.uk/spm/software/)'],...
    cits);
