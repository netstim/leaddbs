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

ea_mkdir(fullfile(options.subj.coregDir, 'transformations'));
movefile(transform{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'brainsfit.mat']);
movefile(transform{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'brainsfit.mat']);

% convert BRAINSFit matrices to 4x4
load([options.subj.coreg.transform.CT.forwardBaseName, 'brainsfit.mat'])
tmat = ea_antsmat2mat(AffineTransform_double_3_3,fixed);
save([options.subj.coreg.transform.CT.forwardBaseName, 'brainsfit44.mat'],'tmat')
load([options.subj.coreg.transform.CT.inverseBaseName, 'brainsfit.mat'])
tmat = ea_antsmat2mat(AffineTransform_double_3_3,fixed);
save([options.subj.coreg.transform.CT.inverseBaseName, 'brainsfit44.mat'],'tmat')

disp('Coregistration done.');
