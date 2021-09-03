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

ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);
ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);
ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);
ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
          [options.subj.postopAnat.(options.subj.postopModality).preproc],...
          [options.subj.postopAnat.(options.subj.postopModality).coreg]);

disp('Coregistration done.');

