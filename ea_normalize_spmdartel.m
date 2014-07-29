function varargout=ea_normalize_spmdartel(options)
% This is a function that normalizes both a copy of transversal and coronar
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
% The procedure used here uses the SPM DARTEL approach to map a patient's
% brain to MNI space directly. Unlike the usual DARTEL-approach, which is
% usually used for group studies, here, DARTEL is used for a pairwise
% co-registration between patient anatomy and MNI template. It has been
% shown that DARTEL also performs superior to many other normalization approaches 
% also  in a pair-wise setting e.g. in 
%   Klein, A., et al. (2009). Evaluation of 14 nonlinear deformation algorithms 
%   applied to human brain MRI registration. NeuroImage, 46(3), 786?802. 
%   doi:10.1016/j.neuroimage.2008.12.037
%
% Since a high resolution is needed for accurate DBS localizations, this
% function applies DARTEL to an output resolution of 0.5 mm isotropic. This
% makes the procedure quite slow.

% The function uses some code snippets written by Ged Ridgway.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn



if ischar(options) % return name of method.
    varargout{1}='SPM DARTEL nonlinear [MR/CT]';
    return
end


segmentresolution=0.5; % resolution of the DARTEL-Warps. Setting this value to larger values will generate the usual DARTEL-Workflow.



usecombined=0; % if set, eauto will try to fuse coronar and transversal images before normalizing them.
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

% check if backup files exist, if not backup

if ~exist([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],'file')
    try    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized]); end
else
    try copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]); end
end
if ~exist([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized],'file')
    try    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized]); end
else
    try    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized]); end
end
if ~exist([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized],'file')
    try    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized]); end
else
    try    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized]); end
end


% First, do the coreg part:

    ea_coreg(options,options.prefs.normalize.coreg);




% now dartel-import the preoperative version.

%try
    disp('Segmenting preoperative version (Import to DARTEL-space)');
    load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob']);
    job.channel.vols{1}=[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized];
    for tpm=1:6
        job.tissue(tpm).tpm=fullfile(fileparts(which('spm')),'toolbox','Seg',['TPM.nii,',num2str(tpm)]);
        if tpm<4
        job.tissue(tpm).native=[0,1];
        else
            job.tissue(tpm).native=[0,0];
        end
    end
    job.resolution=segmentresolution; 
    
    ea_spm_preproc_run(job); % exactly the same as the SPM version ("New Segment" in SPM8) but with an increase in resolution to 0.5 mm iso.
    
    
    disp('*** Segmentation of preoperative MRI worked.');

    
    
    % check if darteltemplate is available, if not generate one
    
    if exist([options.earoot,filesep,'templates',filesep,'dartel',filesep,'dartelmni_6.nii'],'file')
        % There is a DARTEL-Template. Check if it will match:
        Vt=spm_vol([options.earoot,filesep,'templates',filesep,'dartel',filesep,'dartelmni_6.nii']);
        Vp=spm_vol([options.root,filesep,options.patientname,filesep,'rc1pre_tra.nii']);
        if ~isequal(Vp.dim,Vt(1).dim) || ~isequal(Vp.mat,Vt(1).mat) % Dartel template not matching. -> create matching one.
            ea_create_mni_darteltemplate([options.root,filesep,options.patientname,filesep,'rc1pre_tra.nii']);
        end
        
    else % no dartel template present. -> Create matching dartel templates from highres version.
        ea_create_mni_darteltemplate([options.root,filesep,options.patientname,filesep,'rc1pre_tra.nii']);
        
    end
    
    %
    
    


% Normalize to MNI using DARTEL.
matlabbatch{1}.spm.tools.dartel.warp1.images = {
                                                {[options.root,options.prefs.patientdir,filesep,'rc1',options.prefs.prenii_unnormalized,',1']}
                                                {[options.root,options.prefs.patientdir,filesep,'rc2',options.prefs.prenii_unnormalized,',1']}
                                                {[options.root,options.prefs.patientdir,filesep,'rc3',options.prefs.prenii_unnormalized,',1']}
                                                }';
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_5.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_4.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_3.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_2.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_1.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;
jobs{1}=matlabbatch;
%try
    cfg_util('run',jobs);
    disp('*** Dartel coregistration of preoperative version worked.');
%catch
%    error('*** Dartel coregistration failed.');
%end
clear matlabbatch jobs;

% export normalization parameters:
for inverse=0:1
    if inverse
        addstr='_inv';
    else
        addstr='';
    end
    matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
    matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1-inverse 0+inverse];
    matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
    matlabbatch{1}.spm.util.defs.ofname = ['ea',addstr,'_normparams'];
    matlabbatch{1}.spm.util.defs.fnames = '';
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.defs.interp = 1;
    jobs{1}=matlabbatch;

    cfg_util('run',jobs);
    disp('*** Exported normalization parameters to y_ea_normparams.nii');
    clear matlabbatch jobs;
end


% create warped:
% ...pre_tra
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized]}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
try
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of preoperative version worked.');
end
clear matlabbatch jobs;

% ...tra
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
try
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of postoperative version worked.');
end
clear matlabbatch jobs;

% ... cor/sag
try
for fina=1:length(finas)
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{finas{fina}}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
try
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of postoperative version worked.');
end
clear matlabbatch jobs;
end
end

% ...CT
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered]}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
try
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of postoperative CT worked.');
end
clear matlabbatch jobs;



% make normalization "permanent" and include correct bounding box.




for export=1:5
    switch export
        case 2
            outf=options.prefs.prenii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized,',1'];
            gfina=[options.root,options.prefs.patientdir,filesep,options.prefs.gprenii];
        case 2
            outf=options.prefs.tranii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized,',1'];
            gfina=[options.root,options.prefs.patientdir,filesep,options.prefs.gtranii];
        case 3
            outf=options.prefs.cornii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized,',1'];
            gfina=[options.root,options.prefs.patientdir,filesep,options.prefs.gcornii];
        case 4
            outf=options.prefs.sagnii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized,',1'];
            gfina=[options.root,options.prefs.patientdir,filesep,options.prefs.gsagnii];
        case 5
            outf=options.prefs.ctnii;
            fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_coregistered,',1'];
            gfina=[options.root,options.prefs.patientdir,filesep,options.prefs.gctnii];
    end
    
    % save a backup
    
    if exist([options.root,options.prefs.patientdir,filesep,outf],'file')
       try movefile([options.root,options.prefs.patientdir,filesep,outf],[options.root,options.prefs.patientdir,filesep,'ea_backup',date,num2str(now),outf]); end
    end
    
    matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'templates',filesep,'bb.nii,1']
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
        cfg_util('run',jobs);
    end
    clear matlabbatch jobs;
    
    
    % build global versions of files

    try movefile(fina(1:end-2),gfina); end
end







try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gtranii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gcornii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gsagnii]); end


% restore original files (from coregistration)
try copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]); end
try    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized]); end
try    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized]); end




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

