function varargout=ea_normalize_spmdartel_auto(options)
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
% The procedure used here uses the SPM "Segment" routine and
% is probably the most straight-forward way using SPM.
%
% The function uses some code snippets written by Ged Ridgway.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='SPM DARTEL nonlinear (automatic)';
    return
end


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
    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized]);
else
    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]);
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







% first step, coregistration between transversal and coronar/sagittal versions. on full brain


normlog=zeros(4,1); % log success of processing steps. 4 steps: 1. coreg tra and cor, 2. grand mean normalization 3. subcortical normalization 4. subcortical fine normalization that spares the ventricles.


for export=1:2
    for costfun=1:2
        switch export
            case 1
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,',1'];
            case 2
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,',1'];
        end
        
        
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {fina};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfuns{costfun};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        
        jobs{1}=matlabbatch;
        try
            cfg_util('run',jobs);
            normlog(1)=1;
            disp('*** Coregistration between transversal and coronar versions worked.');
            finas{export}=fina; % assign only if worked.
        catch
            disp('*** Coregistration between transversal and coronar versions failed / Using CT Modality.');
            %error('This normalization cannot be performed automatically with eAuto. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
        end
        clear matlabbatch jobs;
        
        matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']
            fina};
        jobs{1}=matlabbatch;
        try % CT
            cfg_util('run',jobs);
            
            
        end
                    clear matlabbatch jobs;

    end
end




% second step, coreg post to pre.

normlog=zeros(4,1); % log success of processing steps. 4 steps: 1. coreg tra and cor, 2. grand mean normalization 3. subcortical normalization 4. subcortical fine normalization that spares the ventricles.



for costfun=1:2
    
    
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']};
    try
        matlabbatch{1}.spm.spatial.coreg.estimate.other = finas;
    catch
        matlabbatch{1}.spm.spatial.coreg.estimate.other={''};
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfuns{costfun};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    jobs{1}=matlabbatch;
    try
        cfg_util('run',jobs);
        normlog(1)=1;
        disp('*** Coregistration between pre and post versions worked.');
    catch
        disp('*** Coregistration between pre and post versions failed.');
        %error('This normalization cannot be performed automatically with LEAD. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
    end
    clear matlabbatch jobs;
    
if exist('finas','var')    
    matlabbatch{1}.spm.util.checkreg.data = [{[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']},finas];
else
    matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']};

end
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;

end



% now dartel-import the preoperative version.

try
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
    
    ea_spm_preproc_run(job,0.5); % exactly the same as the SPM version ("New Segment" in SPM8) but with an increase in resolution to 0.5 mm iso.
    
    
    %cfg_util('run',jobs);
    disp('*** Segmentation of preoperative MRI worked.');
catch
    error('*** Segmentation of preoperative MRI failed.');
end
clear matlabbatch jobs;


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
try
    cfg_util('run',jobs);
    disp('*** Dartel coregistration of preoperative version worked.');
catch
    error('*** Dartel coregistration failed.');
end
clear matlabbatch jobs;



% create warped:
% ...pre_tra
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized]}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of preoperative version worked.');
clear matlabbatch jobs;

% ...tra
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of preoperative version worked.');
clear matlabbatch jobs;

% ... cor/sag
for fina=1:length(finas)
matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{finas{fina}}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
jobs{1}=matlabbatch;    cfg_util('run',jobs);
disp('*** Dartel Warp of preoperative version worked.');
clear matlabbatch jobs;
end




% make normalization "permanent" and include correct bounding box.




for export=2:4
    switch export
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
    end
    
    % save a backup
    
    if exist([options.root,options.prefs.patientdir,filesep,outf],'file')
        movefile([options.root,options.prefs.patientdir,filesep,outf],[options.root,options.prefs.patientdir,filesep,'ea_backup',date,num2str(now),outf]);
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

    movefile(fina(1:end-2),gfina);
end







movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gtranii]);
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gcornii]); end
try movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gsagnii]); end


% restore original files (from coregistration)
copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]);
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

