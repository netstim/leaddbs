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
startdir = ptdir;
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
    elseif strcmp(response,'No')
        disp(['No cortex created in Patient Directory: ' ptdir '/cortex'])
        if ~exist([ptdir 'cortex/CortElecs.mat'],'file')
            qst = {'Do you have subdural electrode coordinates'; 'that you would like to import now?'};
            ImportElecsOption = questdlg(qst,'Import FS'); clear qst
            if strcmp(ImportElecsOption,'Yes')
                vars = {'patientname','ptdir','fsdir'};
                info = load([ptdir '/cortex/CortexHiRes.mat'],vars{:});
                info.uipatdirs = info.ptdir;
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
    elseif strcmp(response,'Yes')
        load([ptdir,'/cortex/CortexHiRes.mat'],'fsdir');
        FsDir = fsdir;
    else
        % Choose Freesurfer Directory
        FsDir = ea_uigetdir(startdir,['Choose Freesurfer Folder for ' patientname]);
    end
end

if ~exist('FsDir','var')
     FsDir = ea_uigetdir(startdir,['Choose Freesurfer Folder for ' patientname]);
end
if iscell(FsDir) && length(FsDir)==1
    FsDir = char(FsDir);
elseif isempty(FsDir)
    disp('No files saved')
    return
elseif ischar(FsDir)
    %nothing
else
    ea_error('Please choose one FS folder at a time')
end

if ~exist([ptdir,'/cortex'],'dir')
    mkdir([ptdir,'/cortex'])
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
AsegFile = [FsDir '/mri/aseg.mgz'];

% Convert T1.mgz to T1.nii (Freesurfer Dependent)
% if ~exist([FsDir '/mri/T1.nii'],'file')
%     system(['mri_convert -i ' FsDir '/mri/T1.mgz ' -o ' FsDir '/mri/T1.nii -it mgz -ot nii'])
% end
    % Notes: need to add PC functionality
    % Notes: need to add ea_libs_helper for Freesurfer compatibility

annot_DKT(1).atlas = 'Desikan-Killiany';
annot_DKT(2).atlas = 'Desikan-Killiany';
annot_DKT(1).filename  = 'label/rh.aparc.annot';
annot_DKT(2).filename  = 'label/lh.aparc.annot';
annot_DKTaseg(1).atlas = 'Desikan-Killiany + Aseg';
annot_DKTaseg(2).atlas = 'Desikan-Killiany + Aseg';
annot_DKTaseg(1).filename  = 'label/rh.aparc.DKTatlas40.annot';
annot_DKTaseg(2).filename  = 'label/lh.aparc.DKTatlas40.annot';
annot_a2009(1).atlas = 'Destrieux';
annot_a2009(2).atlas = 'Destrieux';
annot_a2009(1).filename  = 'label/rh.aparc.a2009s.annot';
annot_a2009(2).filename  = 'label/lh.aparc.a2009s.annot';

% Read Annotation Files
% external/freesurfer/read_annotation.m
for iSide = 1:2
    [annot_DKT(iSide).vert, annot_DKT(iSide).label, annot_DKT(iSide).colortable] = read_annotation([FsDir,filesep,annot_DKT(iSide).filename]);
    [annot_DKTaseg(iSide).vert, annot_DKTaseg(iSide).label, annot_DKTaseg(iSide).colortable] = read_annotation([FsDir,filesep,annot_DKTaseg(iSide).filename]);
    [annot_a2009(iSide).vert, annot_a2009(iSide).label, annot_a2009(iSide).colortable] = read_annotation([FsDir,filesep,annot_a2009(iSide).filename]);
end
annot = annot_DKT;
save([ptdir '/cortex/annot_DKT.mat'],'annot'); clear annot;
annot = annot_DKTaseg;
save([ptdir '/cortex/annot_DKTaseg.mat'],'annot'); clear annot
annot = annot_a2009;
save([ptdir '/cortex/annot_a2009.mat'],'annot'); clear annot


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

%% Translating into the appropriate space
% FS to MRI
CortexHiRes.raw.nativespace = ea_popupquest('In what space was your freesurfer brain reconstructed?',...
    '','Preop','Postop');
disp(['Translating into native ' CortexHiRes.raw.nativespace ' space...'])
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

% Postop to preop
switch CortexHiRes.raw.nativespace
    case 'Postop'
    coregfile = fdir(ptdir,'GenericAffine'); %fdir(ptdir,'2postop');
    % %     ONLY SUPPORTS ANTS COREGISTRATION
    if ~isempty(coregfile) && length(coregfile)==1 && ~isempty(regexp(coregfile,'\wants','match'))
        load([ptdir,filesep,cell2mat(coregfile)])
        try 
            aff = ea_antsmat2mat(AffineTransform_float_3_3,fixed);
        catch
            aff = ea_antsmat2mat(AffineTransform_double_3_3,fixed);
        end
        aff([4,8,12])=aff([13:15]); aff([13:15])=0;
        tform = affine3d(aff);
        CortexHiRes.raw.coregistration = tform;
        CortexHiRes.Vertices_lh = transformPointsInverse(tform, CortexHiRes.Vertices_lh);
        CortexHiRes.Vertices_rh = transformPointsInverse(tform, CortexHiRes.Vertices_rh);
        CortexHiRes.Vertices = transformPointsInverse(tform, CortexHiRes.Vertices);
        
        [CortexHiRes.Vert_mm,CortexHiRes.Vert_vox] = ea_map_coords(CortexHiRes.Vertices','raw_postop_tra.nii',fullfile(ptdir,cell2mat(coregfile)),'anat_t1.nii','ANTS');
    else 
        ea_warning('cannot find corgeistration matrix')
    end
    otherwise
end


%% Save Output to PatientDirectory/cortex/
if ~exist(fullfile(ptdir,'cortex'),'dir')
    mkdir(fullfile(ptdir,'cortex'))
end
    disp(['Saving to ' fullfile(ptdir,'cortex/CortexHiRes.mat') '...'])
    save(fullfile(ptdir,'cortex/CortexHiRes.mat'),'-struct','CortexHiRes')

%     return
    %% Get Hull under dev
%     cmd = '/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/mni/bin:/Applications/freesurfer/tktools';
%     ea_libs_helper(cmd,'PATH')
% 
%     grayfilename = [FsDir '/mri/ribbon.nii'];
%     % if ~exist(grayfilename,'file') && exist([FsDir,'/mri/ribbon.mgz'],'file')
%     %     cmd = sprintf('mri_convert -i %s/mri/ribbon.mgz -o %s/mri/ribbon.nii -it mgz -ot nii',FsDir,FsDir);
%     %     system(cmd)
%     % else ~exist(grayfilename,'file') && ~exist([FsDir,'/mri/ribbon.mgz'],'file')
%     %     ea_warning('Cannot Find mri/ribbon.mgz')
%     % end    
%         %disp('running dbs_gethull.......')
%         %[mask_matrix,mask_indices] = ea_gethull(grayfilename,3,21,.3);
%         save([ptdir '/cortex/hull.mat'],'mask_matrix','mask_indices')
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
    
elseif strcmp(DownsampleOption,'Cancel')
    return

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


function files = fdir(Path,exp)
% Returns contents of path as cell of filenames

if isempty(Path) || ~exist('Path','var')
    Path = pwd;
end
    
contents = dir(Path);
contents = contents(cellfun(@(x) isempty(regexp(x, '^\.', 'once')), {contents.name}));
files = {contents(~[contents.isdir]).name};

if nargin==2
    files = files(~cellfun(@(x) isempty(regexp(x, ['\w' exp],'match')), files));
end
