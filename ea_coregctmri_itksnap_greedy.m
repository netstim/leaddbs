function varargout=ea_coregctmri_itksnap_greedy(options)
% This function uses ITK-SNAPs greedy binary to do an automatic very fast
% rigid registration between post-op CT and MRI.
%
% Some exhaustive search for a good initial configuration is carried out
% before optimsation thus it should perform robust also on heavily tilted
% data (up to 90 degree tilt).
%
% The transformation is stored in ITK .txt format and converted to the
% newer ITK .mat format used by ANTs. Thus all further operations can
% threat it as an ANTs transform.
% __________________________________________________________________________________
% Copyright (C) 2019 University of Luxembourg, Interventional Neuroscience
% Group
% Andreas Husch

if ischar(options) % return name of method (e.g. used in the GUI)
    varargout{1}='ITK-SNAP greedy (automatic)'; %TODO make sure matlab finds the path
    varargout{2}={'SPM8','SPM12'}; % check this
    varargout{3}=['nan'];
    return
end

disp('Running ITK-SNAP greedy rigid registration...');
%% Get itksnap binary (should be in path by following ITKSNAP / Help / Install Comand Line Tools)
if ispc
    greedy = 'greedy.exe';
else
    greedy = '/usr/local/bin/greedy';
end

cmd = [greedy ' '...
   '-d 3 '...
   '-a -dof 6 '...
   '-search 100 90 5 '...
   '-ia-image-centers '...
   '-m MI '...
   '-i ' options.root,options.patientname,filesep,options.prefs.prenii_unnormalized ' ' ...
   options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized ...
   '-o postop_ct2anat_t1_ants1.txt -n 100x50x10'];
   
   

try
    if ~ispc
        system(['bash -c ''' cmd '''']); % running in background (&) is key to avoid crash!
    else
        system(cmd);
    end
  %  while ~exist([options.root,options.patientname,filesep,'postop_ct2anat_t1_ants1.txt'], 'file')
   %     pause(2); % semi busy waiting ... ;-)
   % end
catch
    warning('Could not run ITK-SNAP successfully! Make sure itksnap is available in your path (Launch ITK-SNAP and see in Help / Install Command Line Tools)');
end

%% Save ITK .txt also as ANTS .mat
disp('Generating .mat version of transform to be re-read by LeadDBS for internal use later...');
ea_convertANTsITKTransforms([options.root,options.patientname,filesep,'postop_ct2anat_t1_ants1.txt'], ...
    [options.root,options.patientname,filesep,'postop_ct2anat_t1_ants1.mat']);
ea_convertANTsITKTransforms([options.root,options.patientname,filesep,'postop_ct2anat_t1_ants1.txt'], ...
    [options.root,options.patientname,filesep,'anat_t12postop_ct_ants1.mat'], 1); % generate inverse
disp('Generating .mat done.');

%% Apply registration to generate transformed CT file
disp('Generating image by applying transform...')
ea_apply_coregistration([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized], ...
    [options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized],...
    [options.root,options.patientname,filesep,options.prefs.ctnii_coregistered]);
disp('Coregistration done.');