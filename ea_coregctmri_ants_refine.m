function varargout=ea_coregctmri_ants_refine(options)
% This function uses ANTs to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Advanced Normalization Tools (ANTs) + Subcortical Refine';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');
transform = ea_ants([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.(options.subj.postopModality).preproc],...
    [options.subj.postopAnat.(options.subj.postopModality).coreg],1,{},1,options);

movefile(transform{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'ants.mat']);
movefile(transform{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'ants.mat']);

disp('Coregistration done.');
