function varargout=ea_coregpostopct_brainsfit(options)
% Wrapper function for BRAINSFit registration of post-op CT
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='BRAINSFit (Johnson 2007)';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');

transform = ea_brainsfit([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg]);

movefile(transform{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'brainsfit.mat']);
movefile(transform{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'brainsfit.mat']);

disp('Coregistration done.');
