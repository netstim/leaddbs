function varargout=ea_coregpostopct_ants_multiple(options)
% Wrapper function for ANTs (multiple-run) registration of post-op CT
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='ANTs (Avants 2008) multiple runs';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');
transform1 = ea_ants_linear([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);
transform2 = ea_ants_linear([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);
transform3 = ea_ants_linear([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);
transformFinal = ea_ants_linear([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);

ea_delete([transform1; transform2; transform3]);

movefile(transformFinal{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'ants.mat']);
movefile(transformFinal{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'ants.mat']);

disp('Coregistration done.');
