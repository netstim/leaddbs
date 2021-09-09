function varargout=ea_normalize_spmseg_macaque(options)
% This is a function that normalizes both a copy of transversal and coronal
% images into MNI-space. The goal was to make the procedure both robust and
% automatic, but still, it must be said that normalization results should
% be taken with much care because all reconstruction results heavily depend
% on these results. Normalization of DBS-MR-images is especially
% problematic since usually, the field of view doesn't cover the whole
% brain (to reduce SAR-levels during acquisition) and since electrode
% artifacts can impair the normalization process. Therefore, normalization
% might be best archieved with other tools that have specialized on
% normalization of such image data.
%
% The procedure used here uses the SPM "Segment" routine and
% is probably the most straight-forward way using SPM.
%
% The function uses some code snippets written by Ged Ridgway.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='SPM12 Old Segment [MR/CT]';
    varargout{2}={'SPM12'};
    return
end

warning('off');
usecombined=0; % if set, LEAD will try to fuse coronal and transversal images before normalizing them.
usesegmentnew=0;
costfuns={'nmi','mi','ecc','ncc'};

if exist([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,'.gz'],'file')
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,'.gz']);
    end


    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,'.gz']);
    end
end


% do a linear coregistration into mni space
if options.modality==1 %MR
    expdo=2:4;
elseif options.modality==2 % CT
    expdo=6;
end

% apply estimated transformations to cor and tra.
cnt=1;
for export=expdo
    switch export
        case 2
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized];
                cnt=cnt+1;
            end
        case 3
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized];
                cnt=cnt+1;
            end
        case 4
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized];
                cnt=cnt+1;
            end
        case 6 % CT
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered];
                cnt=cnt+1;
            end
    end
end

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(options.earoot,'toolbox','macaque','templates','mni_hires_t2.nii,1')};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = finas;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',{matlabbatch});


% now segment the preoperative version.
matlabbatch{1}.spm.spatial.preproc.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 0;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
    fullfile(options.earoot,'toolbox','macaque','templates','mni_inia19_prob_c1.nii');
    fullfile(options.earoot,'toolbox','macaque','templates','mni_inia19_prob_c2.nii');
    fullfile(options.earoot,'toolbox','macaque','templates','mni_inia19_prob_c3.nii')
    };
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2; 2; 2; 4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni'; %'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};

jobs{1}=matlabbatch;
try
    spm_jobman('run',jobs);
    disp('*** Segmentation of preoperative version worked.');
catch
    disp('*** Segmentation of transversal version failed.');
    ea_error('This normalization cannot be performed automatically with LEAD. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
end
clear matlabbatch jobs;

if options.modality==1 %MR
    expdo=2:5;
elseif options.modality==2 % CT
    expdo=5:6;
end

% apply estimated transformations to cor and tra.
for export=expdo
    switch export
        case 2
            outf=options.prefs.tranii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized];
        case 3
            outf=options.prefs.cornii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized];
        case 4
            outf=options.prefs.sagnii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized];
        case 5
            outf=options.prefs.prenii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized];
        case 6 % CT
            outf=options.prefs.ctnii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_coregistered];
    end

[~,nm]=fileparts(options.prefs.prenii_unnormalized); % cut off file extension

voxi=[0.22 0.22 0.22]; % export highres
bbi=[-22 20 22; 22 -40 -10]; % with small bounding box

matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[options.root,options.prefs.patientdir,filesep,nm,'_seg_sn.mat']};
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = voxi;
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = bbi;
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fina};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[options.root,options.prefs.patientdir]};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';

jobs{1}=matlabbatch;
try
spm_jobman('run',jobs);
end
clear matlabbatch jobs;
end

% make normalization "permanent" and include correct bounding box.
for export=expdo
    switch export
        case 2
            outf=options.prefs.tranii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized,',1'];
        case 3
            outf=options.prefs.cornii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized,',1'];
        case 4
            outf=options.prefs.sagnii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized,',1'];
        case 5
            outf=options.prefs.prenii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized,',1'];
        case 6 % CT
            outf=options.prefs.ctnii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_coregistered,',1'];
    end


    matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'toolbox',filesep,'macaque',filesep,'templates',filesep,'bb.nii,1'];
        fina
        };
    matlabbatch{1}.spm.util.imcalc.output = outf;
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = ['i2'];
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    jobs{1}=matlabbatch;
    try
        spm_jobman('run',jobs);
    end
    clear matlabbatch jobs;
end

% build global versions of files
for export=expdo
    switch export
        case 2
            outf=options.prefs.gtranii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized];
        case 3
            outf=options.prefs.gcornii;
                        fina=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized];
        case 4
            outf=options.prefs.gsagnii;
                        fina=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized];
        case 5
            outf=options.prefs.gprenii;
            fina=[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized];
        case 6 % CT
            outf=options.prefs.gctnii;
                        fina=[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered];
            wfina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_coregistered];
    end

    [~,nm]=fileparts(options.prefs.prenii_unnormalized); % cut off file extension

    voxi=[0.5 0.5 0.5]; % export highres
    bbi=nan(2,3);

    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[options.root,options.prefs.patientdir,filesep,nm,'_seg_sn.mat']};
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = voxi;
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = bbi;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fina};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[options.root,options.prefs.patientdir]};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';

    jobs{1}=matlabbatch;
    try
        spm_jobman('run',jobs);
        movefile(wfina,outf);
    end
    clear matlabbatch jobs;

end

% export (generalized) normalization parameters:
for inverse=0:1
    if inverse
        addstr='_inv';
    else
        addstr='';
    end

    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[options.root,options.prefs.patientdir,filesep,nm,'_seg',addstr,'_sn.mat']};
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = nan(1,3);
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = nan(2,3);
    matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = ['ea',addstr,'_normparams'];
    matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {[options.root,options.prefs.patientdir,filesep]};
    jobs{1}=matlabbatch;

    spm_jobman('run',jobs);
    disp('*** Exported normalization parameters to y_ea_normparams.nii');
    clear matlabbatch jobs;
end

try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gtranii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gprenii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gsagnii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gcornii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_coregistered],[options.root,options.prefs.patientdir,filesep,options.prefs.gctnii]); end


function [flirtmat spmvoxmat fslvoxmat] = worldmat2flirtmat(worldmat, src, trg)
%worldmat2flirtmat: convert NIfTI world (mm) coordinates matrix to flirt
%
% Example:
%  [flirtmat spmvoxmat fslvoxmat] = worldmat2flirtdmat(worldmat, src, trg);
%
% See also: flirtmat2worldmat, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

if ischar(src)
    src = nifti(src);
end
if ischar(trg)
    trg = nifti(trg);
end

spmvoxmat = inv(src.mat) * worldmat * trg.mat;
addone = eye(4); addone(:, 4) = 1;
fslvoxmat = inv(addone) * spmvoxmat * addone;
trgscl = nifti2scl(trg);
srcscl = nifti2scl(src);
flirtmat = inv( srcscl * fslvoxmat * inv(trgscl) );


function flirtmat_write(fname, mat)
%flirtmat_write: save a 4-by-4 matrix to a file as handled by flirt -init
% Example:
%  flirtmat_write(fname, affinemat)
% See also: flirtmat_read, flirtmat2worldmat, worldmat2flirtmat

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

fid = fopen(fname, 'w');
% Note that transpose is needed to go from MATLAB's column-major matrix to
% writing the file in a row-major way (see also flirtmat_read)
str = sprintf('%f  %f  %f  %f\n', mat');
fprintf(fid, str(1:end-1)); % drop final newline
fclose(fid);


function [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg)
%flirtmat2worldmat: convert saved flirt matrix to NIfTI world coords matrix
% flirt matrix is from text file specified in "flirt -omat mat.txt" command
% world matrix maps from NIfTI world coordinates in target to source. Note:
% mat.txt contains a mapping from source to target in FSL's *scaled* coords
% which are not NIfTI world coordinates, and note src-trg directionality!
% worldmat from this script reproduces "img2imgcoord -mm ...".
%
% The script can also return a matrix to map from target to source voxels
% in MATLAB/SPM's one-based convention, or in FSL's zero-based convention
%
% Example:
%  [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
%
% See also: worldmat2flirtmat, flirtmat_read, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

if ischar(src)
    src = nifti(src);
end
if ischar(trg)
    trg = nifti(trg);
end

if ischar(flirtmat)
    flirtmat = flirtmat_read(flirtmat);
end

% src = inv(flirtmat) * trg
% srcvox = src.mat \ inv(flirtmat) * trg.mat * trgvox
% BUT, flirt doesn't use src.mat, only absolute values of the
% scaling elements from it,
% AND, if images are not radiological, the x-axis is flipped, see:
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0810&L=FSL&P=185638
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0903&L=FSL&P=R93775
trgscl = nifti2scl(trg);
srcscl = nifti2scl(src);
fslvoxmat = inv(srcscl) * inv(flirtmat) * trgscl;

% AND, Flirt's voxels are zero-based, while SPM's are one-based...
addone = eye(4); addone(:, 4) = 1;
spmvoxmat = addone * fslvoxmat * inv(addone);

worldmat = src.mat * spmvoxmat * inv(trg.mat);


%%
function scl = nifti2scl(N)
% not sure if this is always correct with rotations in mat, but seems okay!
scl = diag([sqrt(sum(N.mat(1:3,1:3).^2)) 1]);
if det(N.mat) > 0
    % neurological, x-axis is flipped, such that [3 2 1 0] and [0 1 2 3]
    % have the same *scaled* coordinates:
    xflip = diag([-1 1 1 1]);
    xflip(1, 4) = N.dat.dim(1)-1; % reflect about centre
    scl = scl * xflip;
end
