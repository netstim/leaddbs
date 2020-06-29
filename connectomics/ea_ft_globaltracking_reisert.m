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

addpath(genpath([options.earoot,'ext_libs',filesep,'Fibertools']));
addpath(genpath([options.earoot,'ext_libs',filesep,'mrTools']));
rmpath(genpath([options.earoot,'ext_libs',filesep,'mesoft']));


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
redo=ea_prepare_dti(options);
ea_prepare_hardi(options,redo)

%% mask for tracking
directory=[options.root,options.patientname,filesep];

% 'new segment' options.prefs.prenii_unnormalized
if ~exist([directory,'trackingmask.nii'],'file') || redo
    ea_gentrackingmask(options,0)
end


%% do tracking
matlabbatch{1}.dtijobs.tracking.GTtrack.fname.filenameHARDI = {[directory,options.prefs.HARDI]};
matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.dir = {directory};
matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.fname = 'lowFTR.mat';
matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.masknii.filenameMASKnii = {[directory,'trackingmask.nii']};
matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.masknii.thresholdMask = 0.5;
matlabbatch{1}.dtijobs.tracking.GTtrack.parameters = 1;
matlabbatch{1}.dtijobs.tracking.GTtrack.para_weight.custom_para_weight = para_weight;
matlabbatch{1}.dtijobs.tracking.GTtrack.para_other.custom_para_other = para_other;
matlabbatch{1}.dtijobs.tracking.GTtrack.minlen = 3;
matlabbatch{1}.dtijobs.tracking.GTtrack.maxlen = Inf;

spm_jobman('run',{matlabbatch});
clear matlabbatch;

%% perform sampling to get larger number of tracts.

matlabbatch{1}.dtijobs.tracking.GTtrack_accum.filenameFTR = {[directory,'lowFTR.mat']};
matlabbatch{1}.dtijobs.tracking.GTtrack_accum.fname = options.prefs.FTR_unnormalized;
matlabbatch{1}.dtijobs.tracking.GTtrack_accum.numsamps = 10;
matlabbatch{1}.dtijobs.tracking.GTtrack_accum.numits = 1000000;
matlabbatch{1}.dtijobs.tracking.GTtrack_accum.temp = 0.1;

spm_jobman('run',{matlabbatch});
clear matlabbatch;


% ftr=fiberGT_tool('createEFTR',10^7,10,0.1);
% ftrstruct_write(ftr,[directory,options.prefs.FTR_unnormalized]);

%% export .trk copy for trackvis visualization
ea_ftr2trk([directory,options.prefs.FTR_unnormalized],[directory,options.prefs.b0]); % export unnormalized ftr to .trk
disp('Done.');

%% add methods dump:
cits={
    'Reisert, M., Mader, I., Anastasopoulos, C., Weigel, M., Schnell, S., & Kiselev, V. (2011). Global fiber reconstruction becomes practical. NeuroImage, 54(2), 955?962. http://doi.org/10.1016/j.neuroimage.2010.09.016'
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    };
ea_methods(options,['A whole-brain fiber-set was estimated based using the Gibbs'' tracking approach (Reisert 2011) using standard parameters.',...
    ' This was done within a white-matter mask that was estimated on the anatomical scan using the Unified Segmentation approach (Ashburner 2005) as implemented in ',spm('ver'),' and linearly co-registered to the b0-weighted series.'],cits);
