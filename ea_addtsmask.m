function ea_addtsmask(options)
% creates refining basal ganglia masks in subject directory
if ~exist([options.root,options.patientname,filesep,'bgmsk.nii'],'file')
    disp('Registering subcortical mask (Schoenecker 2008) to subject preoperative anatomy...');
    
    
    copyfile([ea_space(options,'space'),options.primarytemplate,'.nii'],[options.root,options.patientname,filesep,'tmp.nii']);
    copyfile([ea_space(options,'subcortical'),'secondstepmask.nii'],[options.root,options.patientname,filesep,'bgmsk.nii']);
    % copyfile([ea_getearoot,'templates',filesep,'schoenecker',filesep,'thirdstepmask.nii'],[options.root,options.patientname,filesep,'bgmsk2.nii']);
    
    if strfind(options.coregmr.method,' + Subcortical Refine')
        options.coregmr.method=strrep(options.coregmr.method,' + Subcortical Refine','');
    end
    
    ea_coreg2images(options,[options.root,options.patientname,filesep,'tmp.nii'],...
        [options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
        [options.root,options.patientname,filesep,'tmp.nii'],...
        {[options.root,options.patientname,filesep,'bgmsk.nii']},0);
    delete([options.root,options.patientname,filesep,'tmp.nii']);
end
