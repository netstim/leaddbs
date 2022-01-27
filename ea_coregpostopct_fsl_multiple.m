function varargout=ea_coregpostopct_fsl_multiple(options)
% Wrapper function for FSL FLIRT (multiple-run) registration of post-op CT
% __________________________________________________________________________________
% Copyright (C) 2019 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='FLIRT (Jenkinson 2001 & 2002) multiple runs';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');
transform1 = ea_flirt([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);
transformFinal = ea_flirt([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);

ea_delete(transform1);

movefile(transformFinal{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'flirt.mat']);
movefile(transformFinal{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'flirt.mat']);

disp('Coregistration done.');
