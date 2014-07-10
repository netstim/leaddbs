function varargout=ea_normalize_spmnewseg(options)
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
% The procedure used here uses the SPM "New Segment" routine and
% is probably the most straight-forward way using SPM.
%
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='SPM New Segment nonlinear';
    return
end

segmentresolution=0.5; % resolution of the New Segment output.


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



% First, do the coreg part:

ea_coreg(options,options.prefs.normalize.coreg);






% now segment the preoperative version.

disp('Segmenting preoperative version.');
load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob']);
job.channel.vols{1}=[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized];


tpminf=fullfile(fileparts(which('spm')),'toolbox','Seg',['TPM.nii']);
tpmoutf=[options.earoot,'templates',filesep,'TPM.nii'];
if ~exist(tpmoutf,'file')

reslice_nii(tpminf,tpmoutf,[segmentresolution,segmentresolution,segmentresolution],3);

end

for tpm=1:6
    job.tissue(tpm).tpm=[tpmoutf,',',num2str(tpm)];
    if tpm<4
        job.tissue(tpm).native=[0,0];
    else
        job.tissue(tpm).native=[0,0];
        %job.tissue(tpm).warped=[0,1]; % to export c1, c2, c3 images as well.
    end
end
    job.resolution=segmentresolution;
job.warp.write=[1,1]; % export deformation fields.

ea_spm_preproc_run(job); % exactly the same as the SPM version ("New Segment" in SPM8) but with an increase in resolution to 0.5 mm iso.


disp('*** Segmentation of preoperative MRI done.');



% Rename deformation fields:

movefile([options.root,options.patientname,filesep,'y_',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'y_ea_normparams.nii']);
movefile([options.root,options.patientname,filesep,'iy_',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);

% Apply estimated deformation to (coregistered) post-op data.
matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_normparams.nii']};
matlabbatch{1}.spm.util.defs.ofname = '';

postops={options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized,options.prefs.prenii_unnormalized,options.prefs.ctnii_coregistered};
cnt=1;
for postop=1:length(postops)
    if exist([options.root,options.patientname,filesep,postops{postop}],'file')
    matlabbatch{1}.spm.util.defs.fnames{cnt}=[options.root,options.patientname,filesep,postops{postop},',1'];
    cnt=cnt+1;
    end
end

matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.patientname,filesep]};
matlabbatch{1}.spm.util.defs.interp = 1;
jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;


% rename files:

try copyfile([options.root,options.patientname,filesep,'w',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gprenii]); end
try movefile([options.root,options.patientname,filesep,'w',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.prenii]); end
try copyfile([options.root,options.patientname,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gtranii]); end
try movefile([options.root,options.patientname,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii]); end
try copyfile([options.root,options.patientname,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gcornii]); end
try movefile([options.root,options.patientname,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii]); end
try copyfile([options.root,options.patientname,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gsagnii]); end
try movefile([options.root,options.patientname,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii]); end
try copyfile([options.root,options.patientname,filesep,'w',options.prefs.ctnii_coregistered],[options.root,options.patientname,filesep,options.prefs.gctnii]); end
try movefile([options.root,options.patientname,filesep,'w',options.prefs.ctnii_coregistered],[options.root,options.patientname,filesep,options.prefs.ctnii]); end





