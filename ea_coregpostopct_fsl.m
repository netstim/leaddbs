function varargout=ea_coregpostopct_fsl(options)
% Wrapper function for FSL FLIRT registration of post-op CT
% __________________________________________________________________________________
% Copyright (C) 2019 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='FLIRT (Jenkinson 2001 & 2002)';
    return
end

disp('Coregistering post-op CT to pre-op MRI...');
transform = ea_flirt([options.subj.preopAnat.(options.subj.AnchorModality).coreg],...
    [options.subj.postopAnat.CT.preproc],...
    [options.subj.postopAnat.CT.coreg],1);

ea_mkdir(fullfile(options.subj.coregDir, 'transformations'));
movefile(transform{1}, [options.subj.coreg.transform.CT.forwardBaseName, 'flirt.mat']);
movefile(transform{2}, [options.subj.coreg.transform.CT.inverseBaseName, 'flirt.mat']);

% convert affinefile from txt to tmat
tmat = readmatrix([options.subj.coreg.transform.CT.forwardBaseName, 'flirt.mat'],'FileType','text');
save([options.subj.coreg.transform.CT.forwardBaseName, 'flirt44.mat'],'tmat');
tmat = readmatrix([options.subj.coreg.transform.CT.inverseBaseName, 'flirt.mat'],'FileType','text');
save([options.subj.coreg.transform.CT.inverseBaseName, 'flirt44.mat'],'tmat');

%% transform FLIRT coregistration matrix to world matrix CT->native and native->CT
[tmat, ~, ~] = flirtmat2worldmatPaCER(str2num(fileread([options.subj.coreg.transform.CT.forwardBaseName, 'flirt.mat'])),...
    [options.subj.postopAnat.CT.preproc], [options.subj.preopAnat.(options.subj.AnchorModality).coreg], 0);

save([options.subj.coreg.transform.CT.inverseBaseName 'world.mat'],'tmat')
tmat = inv(tmat);
save([options.subj.coreg.transform.CT.forwardBaseName 'world.mat'],'tmat')

disp('Coregistration done.');
