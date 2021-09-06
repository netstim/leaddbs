function varargout=ea_coregctmri_brainsfit_multiple(options)
% This function uses the BRAINSfit to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='BRAINSFit - multiple runs';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');

transform1 = ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);
transform2 = ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);
transform3 = ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);
transformFinal = ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);

ea_delete([transform1; transform2; transform3]);

movefile(transformFinal{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'brainsfit.h5']);
movefile(transformFinal{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'brainsfit.h5']);

disp('Coregistration done.');
