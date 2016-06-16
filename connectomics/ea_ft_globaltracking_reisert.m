function varargout=ea_ft_globaltracking_reisert(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Global Fibertracking (Reisert et al. 2011)';
    varargout{2}={'SPM8','SPM12'};
    return
end

gdti_trackingparams='standard'; % select param-preset here (see below, also to create your own)

try % make sure only one DTI tracker is being used.
    rmpath(genpath([options.earoot,'dev',filesep,'mesoft']));
end

 load([options.root,options.patientname,filesep,options.prefs.bval]);
   [~,bvfname]=fileparts(options.prefs.bval);
   bvals=eval(bvfname);
   bval=max(bvals)/1000;

switch gdti_trackingparams

    case 'hd_book'
        para_weight = 0.0006;
        para_other = [1, 0.001, 50, 5000000000, 0.5, 3.75, 0.2, bval,0.5]';
    case 'hd_a'
        para_weight = 0.006;
        para_other = [1, 0.001, 50, 5000000000, 0.5, 3.75, 0.2, bval,0.5]';
    case 'hd_a_light'
        para_weight = 0.02;
        para_other = [1, 0.001, 50, 300000000, 0.5, 3, 0.2, bval,0.5]';
    case 'hd_a_verylight'
        para_weight = 0.03;
        para_other = [0.1, 0.001, 50, 300000000, 0.5, 3, 0.2, bval,0.5]';
    case 'standard'
        para_weight = 0.058;
        para_other = [0.1, 0.001, 50, 300000000, 1, 3, 0.2, bval,0.5]';
    case 'standard_enhanced'
        para_weight = 0.058;
        para_other = [0.1, 0.001, 50, 300000000, 1, 3, 0.2, bval,0.5]';
        
    case 'sparse'
        para_weight = 0.136;
        para_other = [0.1, 0.001, 50, 5*10^7, 1, 3, 0.2, bval,0.5]';

end



%% build DTD (tensor calculation)


ea_prepare_dti(options)
ea_prepare_hardi(options)

%% mask for tracking

directory=[options.root,options.patientname,filesep];


% 'new segment' options.prefs.prenii_unnormalized
if ~exist([directory,'trackingmask.nii'],'file');
    ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
    
    
    %% Coreg options.prefs.prenii_unnormalized to b0 (for label.mat and FTR-Normalization)
    
    copyfile([directory,options.prefs.prenii_unnormalized],[directory,'c',options.prefs.prenii_unnormalized]);
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,'c',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
        [directory,'c2',options.prefs.prenii_unnormalized,',1']
        };
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [2 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'rb0';
    
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;
    
    
    movefile([directory,'rb0c2',options.prefs.prenii_unnormalized],[directory,'trackingmask.nii']);
    
end
%% do tracking

matlabbatch{1}.dtijobs.tracking.GTtrack.fname.filenameHARDI = {[directory,options.prefs.HARDI]};
matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.dir = {directory};
matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.fname = options.prefs.FTR_unnormalized;
matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.masknii.filenameMASKnii = {[directory,'trackingmask.nii']};
matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.masknii.thresholdMask = 0.5;
matlabbatch{1}.dtijobs.tracking.GTtrack.parameters = 1;
matlabbatch{1}.dtijobs.tracking.GTtrack.para_weight.custom_para_weight = para_weight;
matlabbatch{1}.dtijobs.tracking.GTtrack.para_other.custom_para_other = para_other;
matlabbatch{1}.dtijobs.tracking.GTtrack.minlen = 3;
matlabbatch{1}.dtijobs.tracking.GTtrack.maxlen = Inf;


jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;





%% export .trk copy for trackvis visualization

dnii=ea_load_nii([directory,options.prefs.b0]);
niisize=size(dnii.img); % get dimensions of reference template.
specs.origin=[0,0,0];
specs.dim=niisize;
specs.affine=dnii.mat;

[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
disp('Done.');





