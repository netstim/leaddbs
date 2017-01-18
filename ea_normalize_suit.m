function varargout=ea_normalize_suit(options)
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
% The procedure used here uses the SPM8 "New Segment", or SPM12 "Segment" routine and
% is probably the most straight-forward way using SPM8.
%
% This function uses resize_img.m authored by Ged Rigdway
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='SUIT DARTEL normalization (Diedrichsen 2006)';
    switch spm('ver')
        case 'SPM12'

            if  ~isempty(which('suit_isolate_seg'))
                varargout{2}=1;
            else
                varargout{2}=0;
            end
        otherwise
            varargout{2}=0;
    end
    return
end


directory=[options.root,options.patientname,filesep];

if isfield(options.prefs, 'tranii_unnormalized')
    if exist([directory,options.prefs.tranii_unnormalized,'.gz'],'file')
        try
            gunzip([directory,options.prefs.tranii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.cornii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.sagnii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.prenii_unnormalized,'.gz']);
        end
    end
end

% check for SUIT installation
if isempty(which('suit_isolate_seg')) % this function is only visible while SPM is actually "running" (not just on the path). This needs to happen for SUIT to run.
    warning('Cannot find SUIT, starting SPM12.'); % this should not happen since checked for in the beginning. Still leaving this snippet for robustnes (e.g. if someone closes SPM between starting Lead and pressing run).
    spm fmri
    if isempty(which('suit_isolate_seg')) % still not found.
        ea_error('SUIT toolbox not found. Please install SUIT toolbox for SPM12 first (http://www.diedrichsenlab.org/imaging/suit.htm).');
    end
end

% now segment the preoperative version.
ea_dispt('Isolating cerebellum & brainstem...');




[options,presentfiles]=ea_assignpretra(options);
for an=1:length(presentfiles)
    anats{an}=[directory,presentfiles{an}];
    if exist([directory,'orig_',presentfiles{an}],'file') % restore if SUIT had been run before.
        movefile([directory,'orig_',presentfiles{an}],[directory,presentfiles{an}]);
    end
end




%% 1: isolate cerebellum
suit_defaults;
% temporarily add SPM12s compat folder to path (suit_isolate_seg needs
% spm_load_float).
addpath(fullfile(spm('dir'),'compat'));
suit_isolate_seg(anats);
% restore remove compat from path:
rmpath(fullfile(spm('dir'),'compat'));


%% 2: normalize into MNI using DARTEL
[~,anatbase]=fileparts(presentfiles{1});
matlabbatch{1}.spm.tools.suit.normalise_dartel.subjND.gray = {[directory,anatbase,'_seg1.nii,1']};
matlabbatch{1}.spm.tools.suit.normalise_dartel.subjND.white = {[directory,anatbase,'_seg2.nii,1']};
matlabbatch{1}.spm.tools.suit.normalise_dartel.subjND.isolation = {[directory,'c_',anatbase,'_pcereb.nii,1']};
spm_jobman('run',{matlabbatch});
clear matlabbatch

%% 3: create deformation fields from DARTEL flowfield

% Export normalization parameters:
% backward
switch spm('ver')
    case 'SPM8'
        ea_error('SPM8 is not supported anymore in this version of Lead-DBS');
    case 'SPM12'
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {[directory,'u_a_',anatbase,'_seg1.nii']};
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};
        matlabbatch{1}.spm.util.defs.comp{2}.def = {[ea_space,'suit',filesep,'y_suit2icbm2009b.nii']}; % add deformation from SUIT Dartel space to ICBM 2009b

        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = ['ea_normparams'];
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
end
jobs{1}=matlabbatch;

spm_jobman('run',jobs);
disp('*** Exported normalization parameters to y_ea_normparams.nii');
clear matlabbatch jobs;

% forward (inverse)
switch spm('ver')
    case 'SPM8'
        ea_error('SPM8 is not supported in this version of Lead-DBS anymore');
    case 'SPM12'
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[ea_space,'suit',filesep,'y_icbm2009b2suit.nii']}; % add deformation from ICBM 2009b to SUIT Dartel space
        matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {[directory,'u_a_',anatbase,'_seg1.nii']};
        matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [0 1];
        matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
        matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {''};
        matlabbatch{1}.spm.util.defs.comp{3}.id.space = {[directory,'c_',anatbase,'.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'ea_inv_normparams';
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
end
jobs{1}=matlabbatch;

spm_jobman('run',jobs);
disp('*** Exported normalization parameters to y_ea_inv_normparams.nii');
clear matlabbatch jobs;


ea_apply_normalization(options)
