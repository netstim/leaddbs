function affinefile = ea_docoreg_spm(options,moving,fixed,cfun,doreslice,otherfiles,writeoutmat,refine)
% Wrapper for SPM coregistration
% note: if refine is set, will always reslice.
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
if nargin<6
    refine=0;
end

if refine
    ea_addtsmask(options);
end

% Write out the transform from moving image vox to fixed image mm (also
% from fixed image vox to moving image mm)
if nargin < 6
    writeoutmat = 0;
    affinefile = {''};
else
    if ~writeoutmat
        affinefile={''};
    end
end

% Read the original affine matrix from the moving image (the header will be
% overwritten after spm coregistration.
movingmat = spm_get_space(moving);

fixedmat = spm_get_space(fixed);


if refine % this will apply a subcortical mask after a first whole-brain registration
    % first pass without mask:
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

    spm_jobman('run',{matlabbatch});
    clear matlabbatch
    [pth,movingstem,ext]=fileparts(strrep(moving,',1',''));

    [directory]=fileparts(fixed);
    directory=[directory,filesep];

    % second pass applying subcortical mask:
    disp('Applying subcortical mask (Schoenecker 2008)...');


    matlabbatch{1}.spm.tools.oldnorm.est.subj.source = {[pth,filesep,'r',movingstem,ext]};
    matlabbatch{1}.spm.tools.oldnorm.est.subj.wtsrc = {[directory,'bgmsk.nii,1']};
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.template = {fixed};
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.weight = {[directory,'bgmsk.nii,1']};
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smosrc = 6;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smoref = 6;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.regtype = 'subj';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.cutoff = 25;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits = 0;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.reg = 1;

    spm_jobman('run',{matlabbatch});
    clear matlabbatch


    % apply final transform:
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[pth,filesep,'r',movingstem,'_sn.mat']};
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = [NaN NaN NaN];
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = [NaN NaN NaN
        NaN NaN NaN];

    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[pth,filesep,movingstem,'.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pth};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
    spm_jobman('run',{matlabbatch});
    clear matlabbatch



    % cleanup:
    movefile([pth,filesep,'w',movingstem,ext],[pth,filesep,'r',movingstem,ext]);
    delete([pth,filesep,'r',movingstem,'_sn.mat']);

else

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
        spm_jobman('run',{matlabbatch});
        [pth,movingstem,ext]=fileparts(strrep(moving,',1',''));

    else
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fixed};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {moving};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = otherfiles;
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = cfun;
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [12 10 8 6 4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run',{matlabbatch});

    end


end

[~, mov] = ea_niifileparts(moving);
[~, fix] = ea_niifileparts(fixed);

if writeoutmat
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



%% add methods dump:
cits={
    'Friston, K. J., Ashburner, J. T., Kiebel, S. J., Nichols, T. E., & Penny, W. D. (2011). Statistical Parametric Mapping: The Analysis of Functional Brain Images. Academic Press.'
    };

ea_methods(options,[mov,' was linearly co-registered to ',fix,' using SPM12 (Friston 2011; http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)'],...
    cits);


