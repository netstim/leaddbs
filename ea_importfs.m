function ea_importfs(varargin)
% This function imports freesurfer cortical reconstruction files into 
% patient directory, and converts .pial files to cortex.mat
%
% External Dependencies:
%       mri_convert (Freesurfer)
% 
% Function:
% 1) Locate freesurfer directory
% 2) Save cortex.mat to patientdirectory
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh (UPMC), Brain Modulation Lab
% Ari Kappel

if isfield(varargin{1},'patdir_choosebox')
    
    handles = varargin{1};
    
    if strcmp(handles.patdir_choosebox.String,'Choose Patient Directory')
        ea_error('Please Choose Patient Directory')
    else
        ptdir = handles.patdir_choosebox.String;
    end
    
    tmp = strsplit(handles.patdir_choosebox.String,'/');
    patientname = tmp{end}; clear tmp

elseif isfield(varargin{1},'uipatdirs')
    
    options = varargin{1};
    
    ptdir = options.uipatdirs{1};
    patientname = options.patientname;

end

fsdir = ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname]);
if iscell(fsdir) && length(fsdir)==1
    fsdir = char(fsdir);
else
    ea_error('Please choose one FS folder at a time')
end

% Check fsdir compatibility
if ~exist([fsdir '/mri/T1.mgz'],'file')
    msg = ['Choose Freesurfer Folder for ' patientname ' (Missing ' patientname '_FS/mri/T1.mgz)...'];
    fsdir = ea_uigetdir(ptdir,msg);
elseif ~exist([fsdir '/surf/lh.pial'],'file')
    msg = ['Choose Freesurfer Folder for ' patientname ' (Missing ' patientname '_FS/surf/lh.pial)...'];
    fsdir = ea_uigetdir(ptdir,msg);
elseif~exist([fsdir '/mri/T1.mgz'],'file')
    msg = ['Choose Freesurfer Folder for ' patientname ' (Missing ' patientname '_FS/surf/rh.pial)...'];
    fsdir = ea_uigetdir(ptdir,msg);
end

% Convert T1.mgz to T1.nii (Freesurfer Dependent)
if ~exist([fsdir '/mri/T1.nii'],'file')
    system(['mri_convert -i ' fsdir '/mri/T1.mgz -o ' fsdir '/mri/T1.nii -it mgz -ot nii'])
end
    % Notes: need to add PC functionality
    % Notes: need to add ea_libs_helper for Freesurfer compatibility

cortex.info.patientname = patientname;
cortex.info.ptdir = ptdir;
cortex.info.fsdir = fsdir;
   
disp('Loading reconstruction...')
% Read surface files
[cortex.vert_lh,cortex.tri_lh]= read_surf(fullfile(fsdir,'surf/lh.pial')); % Reading left side pial surface
[cortex.vert_rh,cortex.tri_rh]= read_surf(fullfile(fsdir,'surf/rh.pial')); % Reading right side pial surface

% Generate entire cortex
cortex.vert = [cortex.vert_lh; cortex.vert_rh]; % Combining both hemispheres
cortex.tri = [cortex.tri_lh; (cortex.tri_rh + length(cortex.vert_lh))]; % Combining faces

cortex.tri=cortex.tri+1; % freesurfer starts at 0 for indexing

% Reading in MRI parameters
f=MRIread(fullfile(fsdir,'mri/T1.nii'));

% Translating into the appropriate space
for k=1:size(cortex.vert,1)
    a=f.vox2ras/f.tkrvox2ras*[cortex.vert(k,:) 1]';
    cortex.vert(k,:)=a(1:3)';
end

disp(['Saving to ' fullfile(ptdir,'cortex.mat') '...'])
save(fullfile(ptdir,'cortex.mat'),'cortex')
disp('Done')