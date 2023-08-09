function varargout=ea_coregpostopct_ants(options)
% Wrapper function for ANTs registration of post-op CT
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='ANTs (Avants 2008)';
    return
end

disp('Coregistering postop CT to preop MRI...');
transform = ea_ants_linear(options.subj.preopAnat.(options.subj.AnchorModality).coreg,...
    options.subj.postopAnat.CT.preproc,...
    options.subj.postopAnat.CT.coreg);

ea_mkdir(fullfile(options.subj.coregDir, 'transformations'));
movefile(transform{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'ants.mat']);
movefile(transform{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'ants.mat']);

% convert ANTS matrices to 4x4
load([options.subj.coreg.transform.CT.forwardBaseName, 'ants.mat'])
tmat = ea_antsmat2mat(AffineTransform_float_3_3,fixed);
save([options.subj.coreg.transform.CT.forwardBaseName, 'ants44.mat'],'tmat')
load([options.subj.coreg.transform.CT.inverseBaseName, 'ants.mat'])
tmat = ea_antsmat2mat(AffineTransform_float_3_3,fixed);
save([options.subj.coreg.transform.CT.inverseBaseName, 'ants44.mat'],'tmat')

disp('Coregistration done.');
