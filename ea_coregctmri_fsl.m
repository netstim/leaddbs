function varargout=ea_coregctmri_fsl(options)
% This function uses FSL to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2019 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='FSL FLIRT';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');
transform = ea_flirt([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.(options.subj.postopModality).preproc],...
    [options.subj.postopAnat.(options.subj.postopModality).coreg],1);

movefile(transform{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'flirt.mat']);
movefile(transform{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'flirt.mat']);

disp('Coregistration done.');
