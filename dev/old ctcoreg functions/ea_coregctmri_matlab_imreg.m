function varargout=ea_coregctmri_matlab_imreg(options)
% This function uses the Matlab image toolbox to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Matlab Imreg';
    if exist('edge.m','file') % check for imtbx.
        varargout{2}={'SPM8','SPM12'};
    else
        varargout{2}={};
    end
    varargout{3}=['nan']; % suggestion for alpha-parameter.
    return
end


disp('Loading images...');

% MR
ea_reslice_nii([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'hd_',options.prefs.prenii_unnormalized],[0.5 0.5 0.5],0);
MR=ea_load_nii([options.root,options.patientname,filesep,'hd_',options.prefs.prenii_unnormalized],'simple');

% CT
ea_reslice_nii([options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized],[options.root,options.patientname,filesep,'hd_',options.prefs.rawctnii_unnormalized],[0.5 0.5 0.5],0);
CT=ea_load_nii([options.root,options.patientname,filesep,'hd_',options.prefs.rawctnii_unnormalized],'simple');

disp('Done. Smoothing...');
MR.img(isnan(MR.img))=0;
CT.img(isnan(CT.img))=0;
MR.img(MR.img<100)=0; % remove noise around brain.
sMR=smooth3(MR.img,'gaussian',[15 15 15]);
CT.img(CT.img<0)=0; % remove negative hounsfield parts.
sCT=smooth3(CT.img,'gaussian',[15 15 15]);

disp('Done. Coregistering...');

vizz=1; % visualization on
alphas=options.coregct.coregthreshs;

%% estimate

[optimizer,metric] = imregconfig('multimodal');

RMR  = imref3d(size(MR.img),abs(MR.hdr.mat(6)),abs(MR.hdr.mat(1)),abs(MR.hdr.mat(11)));
RCT = imref3d(size(CT.img),abs(CT.hdr.mat(6)),abs(CT.hdr.mat(1)),abs(CT.hdr.mat(11)));


optimizer.InitialRadius = 0.004;
optimizer.MaximumIterations = 300;
%optimizer.Epsilon=0.000000015;

M = imregtform(sCT,RCT, sMR,RMR, 'affine', optimizer, metric,'PyramidLevels',3,'DisplayOptimization',0);

disp('Done. Writing out rCT image...');
rCT=MR;
rCT.hdr.fname=[options.root,options.patientname,filesep,options.prefs.ctnii_coregistered];
rCT.img = imwarp(CT.img,RCT,M,'bicubic','OutputView',RMR);


spm_write_vol(rCT.hdr,rCT.img);

disp('Done. Cleaning up...');

delete([options.root,options.patientname,filesep,'hd_',options.prefs.prenii_unnormalized]);
delete([options.root,options.patientname,filesep,'hd_',options.prefs.rawctnii_unnormalized]);

disp('Done.');
