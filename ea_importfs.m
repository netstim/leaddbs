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
        ea_error('Please Choose Patient Directory');
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

%% Handle File Options
% Check if CortexHiRes.mat and CortexLowRes_*.mat already exists
files = dir([ptdir '/cortex']); files = files(cellfun(@(x) isempty(regexp(x, '^\.', 'once')), {files.name}));
files = files(~[files.isdir]); files = {files(~cellfun(@isempty , strfind({files.name},'Cortex'))).name};
overwrite = ~cellfun(@isempty,strfind(files,'CortexHiRes.mat'));
overwrite = overwrite+~cellfun(@isempty,strfind(files,'CortexLowRes'));
V = cell(size(overwrite,2)-1);
if size(overwrite,2)>=2
  for f = 1:size(overwrite,2)-1
      tmp = strsplit(files{f+1},'_');
      V{f} = tmp{2}(1:end-4);
  end 
end

if overwrite
    qst = [{'Warning:' }...
        {'There is already a Hi Resolution Cortex defined for this subject.'},...
        {'Are you sure you want to overwrite the previous CortexHiRes.mat?'}];
        %         ['     Found: ' ptdir '/cortex/' filenames{1}]];
    response = questdlg(qst,'Import FreeSurfer folder');
    if strcmp(response,'Cancel')
        disp(['No cortex created in Patient Directory: ' ptdir '/cortex'])
        return
    elseif strcmp(response,'No') || strcmp(response,'Cancel')
        disp(['No cortex created in Patient Directory: ' ptdir '/cortex'])
        if ~exist([ptdir 'cortex/CortElecs.mat'],'file')
            qst = {'Do you have subdural electrode coordinates'; 'that you would like to import now?'};
            ImportElecsOption = questdlg(qst,'Import FS'); clear qst
            if strcmp(ImportElecsOption,'Yes')
                vars = {'patientname','ptdir','fsdir'};
                info = load([ptdir '/cortex/CortexHiRes.mat'],vars{:});
                CortElecs = ea_importcorticalels(info);
                if exist('CortElecs','var')
                    disp(['Saving to ' ptdir '/cortex/CortElecs.mat'])
                    save([ptdir '/cortex/CortElecs.mat'],'-struct','CortElecs')
                end
            elseif strcmp(ImportElecsOption,'No') || strcmp(ImportElecsOption,'Cancel')
                return
            end
        end
        return
    end
end

% Choose Freesurfer Directory
FsDir = ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname]);
if iscell(FsDir) && length(FsDir)==1
    FsDir = char(FsDir);
elseif isempty(FsDir)
    disp('No files saved')
    return
else
    ea_error('Please choose one FS folder at a time')
end

%% Parse Freesurfer Folder

if ~exist([FsDir '/mri/T1.mgz'],'file')
    msg = ['Missing: ' FsDir '/mri/T1.mgz'];
    w = warndlg(['Warning: This may not be a freesurfer folder. ' msg],patientname); waitfor(w);
    FsDir = char(ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname ' (' msg ')']));
elseif ~exist([FsDir '/mri/aseg.mgz'],'file')
    msg = ['Missing: ' FsDir '/mri/aseg.mgz'];
    w = warndlg(['Warning: This may not be a freesurfer folder. ' msg],patientname); waitfor(w);
    FsDir = char(ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname ' (' msg ')']));
elseif ~exist([FsDir '/surf/lh.pial'],'file')
    msg = ['Missing: ' FsDir '/mri/lh.pial'];
    w = warndlg(['Warning: This may not be a freesurfer folder. ' msg],patientname); waitfor(w);
    FsDir = char(ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname ' (' msg ')']));
elseif ~exist([FsDir '/surf/rh.pial'],'file')
    msg = ['Missing: ' FsDir '/mri/rh.pial'];
    w = warndlg(['Warning: This may not be a freesurfer folder. ' msg],patientname); waitfor(w);
    FsDir = char(ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname ' (' msg ')']));
elseif ~exist([FsDir '/label/lh.aparc.a2009s.annot'],'file')
    msg = ['Missing: ' FsDir '/mri/lh.aparc.a2009s.annot'];
    w = warndlg(['Warning: This may not be a freesurfer folder. ' msg],patientname); waitfor(w);
    FsDir = char(ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname ' (' msg ')']));
elseif ~exist([FsDir '/label/rh.aparc.a2009s.annot'],'file')
    msg = ['Missing: ' FsDir '/mri/lh.aparc.a2009s.annot'];
    w = warndlg(['Warning: This may not be a freesurfer folder. ' msg],patientname); waitfor(w);
    FsDir = char(ea_uigetdir(ptdir,['Choose Freesurfer Folder for ' patientname ' (' msg ')']));
end

MriFile  = [FsDir '/mri/T1.mgz'];
LhPial   = [FsDir '/surf/lh.pial'];
RhPial   = [FsDir '/surf/rh.pial'];
% AsegFile = [FsDir '/mri/aseg.mgz'];
% AnnotLH  = [FsDir '/label/lh.aparc.a2009s.annot'];
% AnnotRH  = [FsDir '/label/rh.aparc.a2009s.annot'];
% sAtlas.Name = 'Destrieux';

% Convert T1.mgz to T1.nii (Freesurfer Dependent)
% if ~exist([FsDir '/mri/T1.nii'],'file')
%     system(['mri_convert -i ' FsDir '/mri/T1.mgz ' -o ' FsDir '/mri/T1.nii -it mgz -ot nii'])
% end
    % Notes: need to add PC functionality
    % Notes: need to add ea_libs_helper for Freesurfer compatibility

% Read Annotation Files
% external/freesurfer/read_annotation.m
%  [vertices.lh, label.lh, colortable.lh] = read_annotation(AnnotLH);
%  [vertices.rh, label.rh, colortable.rh] = read_annotation(AnnotRH);
    
%     AnnotLhFiles = {file_find(FsDir, 'lh.pRF.annot', 2), file_find(FsDir, 'lh.aparc.a2009s.annot', 2), file_find(FsDir, 'lh.aparc.annot', 2), file_find(FsDir, 'lh.BA.annot', 2), file_find(FsDir, 'lh.BA.thresh.annot', 2), file_find(FsDir, 'lh.aparc.DKTatlas40.annot', 2), ...
%                 file_find(FsDir, 'lh.PALS_B12_Brodmann.annot', 2), file_find(FsDir, 'lh.PALS_B12_Lobes.annot', 2), file_find(FsDir, 'lh.PALS_B12_OrbitoFrontal.annot', 2), file_find(FsDir, 'lh.PALS_B12_Visuotopic.annot', 2), file_find(FsDir, 'lh.Yeo2011_7Networks_N1000.annot', 2), file_find(FsDir, 'lh.Yeo2011_17Networks_N1000.annot', 2)};
% AnnotRhFiles = {file_find(FsDir, 'rh.pRF.annot', 2), file_find(FsDir, 'rh.aparc.a2009s.annot', 2), file_find(FsDir, 'rh.aparc.annot', 2), file_find(FsDir, 'rh.BA.annot', 2), file_find(FsDir, 'rh.BA.thresh.annot', 2), file_find(FsDir, 'rh.aparc.DKTatlas40.annot', 2), ...
%                 file_find(FsDir, 'rh.PALS_B12_Brodmann.annot', 2), file_find(FsDir, 'rh.PALS_B12_Lobes.annot', 2), file_find(FsDir, 'rh.PALS_B12_OrbitoFrontal.annot', 2), file_find(FsDir, 'rh.PALS_B12_Visuotopic.annot', 2), file_find(FsDir, 'rh.Yeo2011_7Networks_N1000.annot', 2), file_find(FsDir, 'rh.Yeo2011_17Networks_N1000.annot', 2)};
% AnnotLhFiles(cellfun(@isempty, AnnotLhFiles)) = [];
% AnnotRhFiles(cellfun(@isempty, AnnotRhFiles)) = [];

%% Create Hi Resolution Cortex
CortexHiRes.patientname = patientname;
CortexHiRes.ptdir = ptdir;
CortexHiRes.fsdir = FsDir;
   
disp('Loading reconstruction...')
% Read surface files
[CortexHiRes.raw.Vertices_lh,CortexHiRes.raw.Faces_lh]= read_surf(LhPial); % Reading left side pial surface
[CortexHiRes.raw.Vertices_rh,CortexHiRes.raw.Faces_rh]= read_surf(RhPial); % Reading right side pial surface

% Generate entire cortex
CortexHiRes.raw.Vertices = [CortexHiRes.raw.Vertices_lh; CortexHiRes.raw.Vertices_rh]; % Combining both hemispheres
CortexHiRes.raw.Faces = [CortexHiRes.raw.Faces_lh; (CortexHiRes.raw.Faces_rh + length(CortexHiRes.raw.Vertices_lh))]; % Combining Faces

% Reading in MRI parameters
T1nii=MRIread(MriFile);

% Translating into the appropriate space
disp('Translating into native space...')
aff = T1nii.vox2ras/T1nii.tkrvox2ras;
aff([4,8,12])=aff([13:15]); aff([13:15])=0;
tform = affine3d(aff); CortexHiRes.raw.tform = tform;
CortexHiRes.Vertices_lh = transformPointsForward(tform,CortexHiRes.raw.Vertices_lh);
CortexHiRes.Vertices_rh = transformPointsForward(tform,CortexHiRes.raw.Vertices_rh);
CortexHiRes.Vertices = transformPointsForward(tform,CortexHiRes.raw.Vertices);

% freesurfer starts at 0 for indexing
CortexHiRes.Faces_lh=CortexHiRes.raw.Faces_lh+1; 
CortexHiRes.Faces_rh=CortexHiRes.raw.Faces_rh+1;
CortexHiRes.Faces=CortexHiRes.raw.Faces+1;

%% Save Output to PatientDirectory/cortex/
if ~exist(fullfile(ptdir,'cortex'),'dir')
    mkdir(fullfile(ptdir,'cortex'))
end
    disp(['Saving to ' fullfile(ptdir,'cortex/CortexHiRes.mat') '...'])
    save(fullfile(ptdir,'cortex/CortexHiRes.mat'),'-struct','CortexHiRes')

%% Option to Downsample CortexHiRes
% newNbVertices = '15000';
qst = {'Would you like to downsample the high '; sprintf('resolution cortex with %d vertices?',size(CortexHiRes.Vertices,1))};
DownsampleOption = questdlg(qst,'Import FreeSurfer');


if strcmp(DownsampleOption,'Yes')
    
    newNbVertices = inputdlg({'Enter the number of vertices for low resolution cortex surface:'},...
        'Import FreeSurfer folder',1,{'15000'});
    
    if ~isempty(V)
        overwrite = strcmp(V,strcat(newNbVertices,'V'));
        qst = [{'Warning:' },...
            strcat({'There is already a Low Resolution Cortex ('},...
            newNbVertices,{'V) for this subject.'}),...
            {['Are you sure you want to overwrite the previous CortexLowRes_' V{overwrite} '.mat?']}];
        response = questdlg(qst,'Import FreeSurfer folder');
    else 
        response = [];
    end
    if isempty(V) || ~isempty(response) && ~strcmp(response,{'Cancel'}) && ~strcmp(response,{'No'})
    newNbVertices = str2double(newNbVertices);
    oldNbVertices = size(CortexHiRes.Vertices,1);
    
    if isempty(newNbVertices) || isnan(newNbVertices) || newNbVertices==0
        fprintf('Cortex not resampled, Hi Resolution only %d Vertices',size(CortexHiRes.Vertices,1))
    else
        
        if (newNbVertices >= oldNbVertices)
            sprintf('Cortex> Surface has %d vertices, cannot downsample to %d vertices.', oldNbVertices, newNbVertices);
            return;
        end
        
        nVertHemi = newNbVertices/2;
        CortexLowRes.patientname = CortexHiRes.patientname;
        CortexLowRes.ptdir = CortexHiRes.ptdir;
        CortexLowRes.fsdir = CortexHiRes.fsdir;
        
        fprintf('Downsampling Cortex From %d Vertices to %d Vertices...\n',oldNbVertices,newNbVertices)
        [CortexLowRes.Vertices_lh, CortexLowRes.Faces_lh] = ea_downsamplecortex(CortexHiRes.Vertices_lh, CortexHiRes.Faces_lh, nVertHemi, 'reducepath');
        [CortexLowRes.Vertices_rh, CortexLowRes.Faces_rh] = ea_downsamplecortex(CortexHiRes.Vertices_rh, CortexHiRes.Faces_rh, nVertHemi, 'reducepath');
        [CortexLowRes.Vertices, CortexLowRes.Faces] = ea_downsamplecortex(CortexHiRes.Vertices, CortexHiRes.Faces, newNbVertices, 'reducepath');
        
    end
      
    disp(['Saving to ' fullfile(ptdir,['cortex/CortexLowRes_' num2str(newNbVertices) 'V.mat']) '...'])
    save(fullfile(ptdir,['cortex/CortexLowRes_' num2str(newNbVertices) 'V.mat']),'-struct','CortexLowRes')
    end
end

%% Import Cortical Electrodes
% Guarantee Options
options.patientname = patientname;
options.uipatdirs = ptdir;
options.fsdir = FsDir;

qst = {'Do you have subdural electrode coordinates'; 'that you would like to import now?'};
ImportElecsOption = questdlg(qst,'Import FS'); clear qst
if strcmp(ImportElecsOption,'Yes')
    CortElecs = ea_importcorticalels(options);
end

% Save To PatientDirectory/cortex
if exist('CortElecs','var')
disp(['Saving to ' ptdir '/cortex/CortElecs.mat'])
    save([ptdir '/cortex/CortElecs.mat'],'-struct','CortElecs')
end

%%
disp('Done')