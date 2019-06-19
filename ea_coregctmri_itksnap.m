function varargout=ea_coregctmri_itksnap(options)
% This function uses ITK-SNAP to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2019 University of Luxembourg, Interventional Neuroscience
% Group
% Andreas Husch

if ischar(options) % return name of method.
    varargout{1}='ITK-SNAP (ANTs compatible ITK transform file)';
    varargout{2}={'SPM8','SPM12'};
    varargout{3}=['nan']; % suggestion for alpha-parameter.
    return
end

disp('Opening ITK-SNAP for coregistration of postop CT to preop MRI...')
disp('Make sure to save transformation file as "postop_ct2anat_t1_ants1.txt"! Process will continue afer transformation file exists.');
try
    system(['/usr/local/bin/itksnap ' ...
        '-g ' options.root,options.patientname,filesep,options.prefs.prenii_unnormalized ' ' ...
        '-o ' options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized]);
    while ~exist([options.root,options.patientname,filesep,'postop_ct2anat_t1_ants1.txt'], 'file')
        pause(2);
    end
catch
    warning('Coulnt run ITK-SNAP successfully! Make sure itksnap is available in your path (Launch ITK-SNAP and see in Help / Install Command Line Tools)');
end
%% TODO apply regisrtation to generate file: [options.root,options.patientname,filesep,options.prefs.ctnii_coregistered
disp('Generating image by appling transform...')
ea_apply_coregistration([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized], ...
    [options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized],...
    [options.root,options.patientname,filesep,options.prefs.ctnii_coregistered]);
disp('Coregistration done.');