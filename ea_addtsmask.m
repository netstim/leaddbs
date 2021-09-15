function ea_addtsmask(options,both)
% creates refining basal ganglia masks in subject directory
if ~exist('both','var')
    both=0;
end
if ~exist([options.root,options.patientname,filesep,'bgmsk.nii'],'file') || both
    disp('Registering subcortical mask (Schoenecker 2008) to subject preoperative anatomy...');

    copyfile([ea_space(options,'subcortical'),'secondstepmask.nii'],[options.root,options.patientname,filesep,'bgmsk.nii']);

    if both
        copyfile([ea_space(options,'subcortical'),'thirdstepmask.nii'],[options.root,options.patientname,filesep,'bgmsk2.nii']);
    end

    if strfind(options.coregmr.method,' + Subcortical Refine')
        options.coregmr.method=strrep(options.coregmr.method,' + Subcortical Refine','');
    end

    for bt=1:1+both
        if bt==2
            btts='2';
        else
            btts='';
        end
        if exist([options.root,options.patientname,filesep,'glanat.nii'],'file') % if normalization has been done already, use inverse warp to map bgmask to anat_t1.nii
            ea_apply_normalization_tofile(options,{[options.root,options.patientname,filesep,'bgmsk',btts,'.nii']},{[options.root,options.patientname,filesep,'bgmsk',btts,'.nii']},1,0,[options.root,options.patientname,filesep,options.prefs.prenii]);
        else % if not, use a simple linear coregistration
            ea_coregimages(options,[ea_space(options,'space'),options.primarytemplate,'.nii'],...
                [options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
                [options.root,options.patientname,filesep,'tmp.nii'],...
                {[options.root,options.patientname,filesep,'bgmsk.nii']},0);
            ea_delete([options.root,options.patientname,filesep,'tmp.nii']);
        end
    end
end
