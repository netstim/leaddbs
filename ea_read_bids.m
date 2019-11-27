function options = ea_read_bids(options,bids_subject_folder,run)

%% get prefs
if ~exist('options','var') || isempty(options)
    options.prefs = ea_prefs('');
end

%% Choose the BIDS subject folder to be imported and read folder structure
if ~exist('bids_subject_folder','var')
    bids_subject_folder = uigetdir;
end

if strcmp(bids_subject_folder(end),filesep)
    bids_subject_folder(end)=[];
end

[bids_folder,id]=fileparts(bids_subject_folder);

if ~strcmp(id(1:3),'sub')
    error('Not a valid BIDS subject folder')
end

%% Check for sessions 

if ~isfield(options.prefs,'bids_session_preop') 
    options.prefs.bids_session_preop = 'ses-preDBS';
end
if ~isfield(options.prefs,'bids_session_postop') 
    options.prefs.bids_session_postop = 'ses-postDBS';% ADD THIS TO LEAD PREFS FILE
end

preop_folder = dir([bids_subject_folder filesep options.prefs.bids_session_preop '*']);
if ~isempty(preop_folder)
    preop_folder = strcat([bids_subject_folder filesep],{preop_folder([preop_folder(:).isdir]).name}');
else 
    preop_folder = bids_subject_folder;
end


postop_folder = dir([bids_subject_folder filesep options.prefs.bids_session_postop '*']);
if ~isempty(postop_folder)
    postop_folder = strcat([bids_subject_folder filesep],{postop_folder([postop_folder(:).isdir]).name}');
else 
    postop_folder = bids_subject_folder;
end



%% Create a lead derivatives folder
derivatives_folder = fullfile(bids_folder,'derivatives','leaddbs',id);
if ~exist(derivatives_folder,'dir')
    mkdir(derivatives_folder)
end

%% set patientdir to the newly created derivatives folder

options.uipatdirs = {derivatives_folder};

options.prefs.prenii_unnormalized_t2star = 'anat_t2star.nii';
options.prefs.prenii_unnormalized_swi = 'anat_swi.nii';
options.prefs.prenii_unnormalized_fgatir = 'anat_fgatir.nii';

lead2bids_lookup = ...
    {'prenii_unnormalized'          preop_folder    ['anat' filesep '*T2w.nii.gz']
    'prenii_unnormalized_t1'        preop_folder    ['anat' filesep '*T1w.nii.gz']
    'prenii_unnormalized_pd'        preop_folder    ['anat' filesep '*PD.nii.gz']
    'prenii_unnormalized_t2star'    preop_folder    ['anat' filesep '*T2star.nii.gz']
    'prenii_unnormalized_swi'       preop_folder    ['anat' filesep '*SWI.nii.gz']
    'prenii_unnormalized_fgatir'    preop_folder    ['anat' filesep '*FGATIR.nii.gz']
    'tranii_unnormalized'           postop_folder   ['anat' filesep '*tra*T2w.nii.gz'] % needs refinement
    'sagnii_unnormalized'           postop_folder   ['anat' filesep '*sag*T2w.nii.gz'] % needs refinement
    'cornii_unnormalized'           postop_folder   ['anat' filesep '*cor*T2w.nii.gz'] % needs refinement
    'rawctnii_unnormalized'         postop_folder   ['ct' filesep '*ct.nii.gz']
    'rest'                          preop_folder    ['func' filesep '*bold.nii.gz']
    'b0'                            preop_folder    ['dwi' filesep '*b0.nii.gz']
    'fa'                            preop_folder    ['dwi' filesep '*fa.nii.gz']
    'fa2anat'                       preop_folder    ['dwi' filesep '*fa2anat.nii.gz']
    'dti'                           preop_folder    ['dwi' filesep '*dwi.nii.gz']};

for a =1 :size(lead2bids_lookup,1)
    for b = 1:length(lead2bids_lookup{a,2})
        files=[];
        files = dir(fullfile(lead2bids_lookup{a,2}{b},lead2bids_lookup{a,3}));
        
        for c=1:length(files)
            if c==1
                outname = fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.nii.gz']);
            else
                outname = fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '_' num2str(c-1) '.nii.gz']);
                warning(['More than 1 file for ' lead2bids_lookup{a,1} ' found!'])
            end
            
             copyfile(fullfile(files(c).folder,files(c).name),fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.nii.gz']));
             gunzip(outname)
             delete(outname)
    
            if strcmp(lead2bids_lookup{a,1},'dti')
                if exist(fullfile(files(c).folder,[files(c).name(1:end-7) '.bval']),'file')
                    movefile(fullfile(files(c).folder,[files(c).name(1:end-7) '.bval']),fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.bval']));
                else
                    warning('No .bval file found for dMRI image.');
                end
                if exist(fullfile(files(c).folder,[files(c).name(1:end-7) '.bvec']),'file')
                    movefile(fullfile(files(c).folder,[files(c).name(1:end-7) '.bvec']),fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.bvec']));
                else
                    warning('No .bvec file found for dMRI image.');
                end
            end
        end
    end
end

cd(options.uipatdirs{1})

if exist('run','var') && run
%     options.normalize.do = true;
%     options.normalize.check = false;
%     options.coregmr.do = 1;
%     options.coregct.do = true;
%     options.scrf.do = 1;


options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.dicomimp.do = 0;
options.dicomimp.method = 1;
options.assignnii = 0;
options.normalize.do = true;
options.normalize.settings = [];
options.normalize.method = 'ea_normalize_ants';
options.normalize.methodn = 9;
options.normalize.check = false;
options.normalize.refine = 0;
options.coregmr.check = 0;
options.coregmr.method = 'SPM';
options.coregmr.do = 1;
options.overwriteapproved = 0;
options.coregct.do = true;
options.coregct.method = 'ea_coregctmri_ants';
options.coregct.methodn = 7;
options.modality = 1;
options.verbose = 3;
options.sides = [1 2];
options.doreconstruction = false;
options.maskwindow = 10;
options.automask = 1;
options.autoimprove = 0;
options.axiscontrast = 8;
options.zresolution = 10;
options.atl.genpt = false;
options.atl.normalize = 0;
options.atl.can = true;
options.atl.pt = 0;
options.atl.ptnative = false;
options.native = 0;
options.d2.col_overlay = 1;
options.d2.con_overlay = 1;
options.d2.con_color = [1 1 1];
options.d2.lab_overlay = 0;
options.d2.bbsize = 50;
options.d2.backdrop = 'MNI_ICBM_2009b_NLIN_ASYM T1';
options.d2.fid_overlay = 1;
options.d2.write = false;
options.d2.atlasopacity = 0.15;
options.d2.writeatlases = 1;
options.manualheightcorrection = false;
options.scrf.do = 1;
options.scrf.mask = 2;
options.d3.write = false;
options.d3.prolong_electrode = 2;
options.d3.verbose = 'on';
options.d3.elrendering = 1;
options.d3.exportBB = 0;
options.d3.hlactivecontacts = 0;
options.d3.showactivecontacts = 1;
options.d3.showpassivecontacts = 1;
options.d3.showisovolume = 0;
options.d3.isovscloud = 0;
options.d3.mirrorsides = 0;
options.d3.autoserver = 0;
options.d3.expdf = 0;
options.d3.writeatlases = 1;
options.numcontacts = 4;
options.entrypoint = 'STN, GPi or ViM';
options.entrypointn = 1;
options.writeoutpm = 1;
options.elmodeln = 1;
options.elmodel = 'Medtronic 3389';
options.atlasset = 'DISTAL Minimal (Ewert 2017)';
options.atlassetn = 13;
options.reconmethod = 'Refined TRAC/CORE';
options.expstatvat.do = 0;
options.fiberthresh = 10;
options.writeoutstats = 1;
options.leadprod = 'dbs';
ea_run('run',options)
end
