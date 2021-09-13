function varargout=ea_ft_mesotracking_reisert(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Mesoscopic Fibertracking (Konopleva et al. 2018)';
    varargout{2}={'SPM8','SPM12'};
    return
end

gdti_trackingparams='standard'; % select param-preset here (see below, also to create your own)
% make sure only one DTI tracker is being used.
rmpath(genpath([options.earoot,'ext_libs',filesep,'Fibertools']));
rmpath(genpath([options.earoot,'ext_libs',filesep,'mrTools']));
addpath(genpath([options.earoot,'ext_libs',filesep,'mesoft']));

switch gdti_trackingparams

    case 'hd_book'
        para_weight = 0.0006;
        para_other = [1
            0.001
            50
            5000000000
            0.5
            3.75
            0.2
            1];
    case 'hd_a'
        para_weight = 0.006;
        para_other = [1
            0.001
            50
            5000000000
            0.5
            3.75
            0.2
            1];
    case 'hd_a_light'
        para_weight = 0.02;
        para_other = [1
            0.001
            50
            300000000
            0.5
            3
            0.2
            1];
    case 'hd_a_verylight'
        para_weight = 0.03;
        para_other = [0.1
            0.001
            50
            300000000
            0.5
            3
            0.2
            1];
    case 'standard'
        para_weight = 0.058;
        para_other = [0.1
            0.001
            50
            300000000
            1
            3
            0.2
            1];
    case 'standard_enhanced'
        para_weight = 0.058;
        para_other = [0.1
            0.001
            50
            300000000
            1
            3
            0.2
            1.5];

end

directory=[options.root,options.patientname,filesep];
ea_prepare_dti(options)

% create c2 from anat
if ~exist([directory,'trackingmask.nii'],'file');
    ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0);

    %% coreg anat to b0
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,'c2',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',{matlabbatch}); clear matlabbatch

    movefile([directory,'rc2',options.prefs.prenii_unnormalized],[directory,'trackingmask.nii']);
end

% create c1 from anat
if ~exist([directory,'gmmask.nii'],'file');
    ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0);
    %% coreg anat to b0
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,'c1',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',{matlabbatch}); clear matlabbatch

    movefile([directory,'rc1',options.prefs.prenii_unnormalized],[directory,'gmmask.nii']);
end


%% mesoft part goes here
[~,dfn]=fileparts(options.prefs.dti);
ea_delete([directory,dfn,'_FTR.mat']);


init_mesoft;
mesoGT_tool('loadData','nii',[directory,options.prefs.dti],{[directory,options.prefs.bvec],[directory,options.prefs.bval],...
     },{[directory,'gmmask.nii'],[directory,'trackingmask.nii']},[128,128]);
mesoGT_tool('reset');
mesoGT_tool('start');

movefile([directory,dfn,'_FTR.mat'],[directory,options.prefs.FTR_unnormalized]);
delete(findobj('tag','fiberGT_main'))

%% export .trk copy for trackvis visualization
ea_ftr2trk([directory,options.prefs.FTR_unnormalized],[directory,options.prefs.b0]); % export unnormalized ftr to .trk
disp('Done.');

%% add methods dump:
cits={'Konopleva, L., Ilyasov, K. A., Skibbe, H., Kiselev, V. G., Kellner, E., Dhital, B., & Reisert, M. (2018). Modelfree global tractography. NeuroImage. http://doi.org/10.1016/j.neuroimage.2018.03.058'
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    };
ea_methods(options,['A whole-brain fiber-set was estimated based using a model-free implementation of the Meso-tracking approach (Konopleva et al. 2018) using standard parameters.',...
    ' This was done within a white-matter mask that was estimated on the anatomical scan using the Unified Segmentation approach (Ashburner 2005) as implemented in ',spm('ver'),' and linearly co-registered to the b0-weighted series.'],cits);
