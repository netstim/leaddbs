function options = ea_read_bids(options,bids_subject_folder,run)

%% get prefs
if ~exist('options','var') || isempty(options)
    options.prefs = ea_prefs('');
end

%% Choose the BIDS subject folder to be imported and read folder structure
if ~exist('bids_subject_folder','var')
    bids_subject_folder = uigetdir;
end

% remove filesep if required
if strcmp(bids_subject_folder(end),filesep)
    bids_subject_folder(end)=[];
end

% separate bids root folder from id folder and get id:
[bids_folder,id]=fileparts(bids_subject_folder);

% at this point lead checks that the subject folder begins with 'sub', not sure this is necessary:
if ~startsWith(id,'sub')
    error('Not a valid BIDS subject folder')
end

%% Check for sessions -> Potentially add further sessions e.g. previous preop without contrast agent or ct. 

if ~isfield(options.prefs,'bids_session_preop') 
    options.prefs.bids_session_preop = 'ses-preDBS';
end

if ~isfield(options.prefs,'bids_session_postop') 
    options.prefs.bids_session_postop = 'ses-postDBS';% ADD THIS TO LEAD PREFS FILE
end

%% get all preop folders
preop_folder = dir([bids_subject_folder filesep options.prefs.bids_session_preop '*']);
if ~isempty(preop_folder)
    preop_folder = strcat([bids_subject_folder filesep],{preop_folder([preop_folder(:).isdir]).name}');
else 
    preop_folder = {bids_subject_folder};
end

%% Get all postop folders
postop_folder = dir([bids_subject_folder filesep options.prefs.bids_session_postop '*']);
if ~isempty(postop_folder)
    postop_folder = strcat([bids_subject_folder filesep],{postop_folder([postop_folder(:).isdir]).name}');
else 
    postop_folder = {bids_subject_folder};
end



%% Create a lead derivatives folder in the BIDS root
derivatives_folder = fullfile(bids_folder,'derivatives','leaddbs',id);
if ~exist(derivatives_folder,'dir')
    mkdir(derivatives_folder)
end

%% set patientdir to the newly created derivatives folder

options.uipatdirs = {derivatives_folder};

% add the definitions of new sequences to the options for later lookup, add your gadolinium sequence as a new definition
% e.g. options.prefs.prenii_unnormalized_t1gd = 'anat_t1gd.nii';
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
    %here you would add 'prenii_unnormalized_t1gd'      preop_folder    ['anat' filesep '*T1gd.nii.gz']
    'tranii_unnormalized'           postop_folder   ['anat' filesep '*tra*T2w.nii.gz'] % needs refinement
    'sagnii_unnormalized'           postop_folder   ['anat' filesep '*sag*T2w.nii.gz'] % needs refinement
    'cornii_unnormalized'           postop_folder   ['anat' filesep '*cor*T2w.nii.gz'] % needs refinement
    'rawctnii_unnormalized'         postop_folder   ['ct' filesep '*ct.nii.gz']
    'rest'                          preop_folder    ['func' filesep '*bold.nii.gz']
    'b0'                            preop_folder    ['dwi' filesep '*b0.nii.gz']
    'fa'                            preop_folder    ['dwi' filesep '*fa.nii.gz']
    'fa2anat'                       preop_folder    ['dwi' filesep '*fa2anat.nii.gz']
    'dti'                           preop_folder    ['dwi' filesep '*diff.nii.gz']};

for a =1 :size(lead2bids_lookup,1) % run through all lookup entries
    for b = 1:length(lead2bids_lookup{a,2}) % run through the different pre- and postop folders folders (more than one possible)
        files=[];
        % find relevant files based on the look up table here it would be possible to implement a regex match
        files = dir(fullfile(lead2bids_lookup{a,2}{b},lead2bids_lookup{a,3})); % lead2bids_lookup{a,2}{b} = folder, lead2bids_lookup{a,3} = filename
        
        for c=1:length(files) % run through the found files
            if c==1
                outname = fullfile(derivatives_folder,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.nii.gz']); % rename the file for lead dbs
            else
                outname = fullfile(derivatives_folder,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '_' num2str(c-1) '.nii.gz']); % add a running number to the remaining files
                warning(['More than 1 file for ' lead2bids_lookup{a,1} ' found!']) % warn that more than one file is found
            end
            
             copyfile(fullfile(files(c).folder,files(c).name),outname); % copy the file
             
             gunzip(outname,derivatives_folder) % unzip the file 
             delete(outname) % delete the zipped file
    
            if strcmp(lead2bids_lookup{a,1},'dti') % check if it is a dti file and copy the bval and bvec files
                if exist(fullfile(files(c).folder,[files(c).name(1:end-7) '.bval']),'file') % check for presence of bval file
                    copyfile(fullfile(files(c).folder,[files(c).name(1:end-7) '.bval']),fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.bval']));
                else
                    warning('No .bval file found for dMRI image.'); % warn if no bval file is found
                end
                if exist(fullfile(files(c).folder,[files(c).name(1:end-7) '.bvec']),'file') % check if bvec file is present
                    copyfile(fullfile(files(c).folder,[files(c).name(1:end-7) '.bvec']),fullfile(options.prefs.patientdir,[options.prefs.(lead2bids_lookup{a,1})(1:end-4) '.bvec']));
                else
                    warning('No .bvec file found for dMRI image.'); % warn if no bvec file is found
                end
            end
        end
    end
end

cd(options.uipatdirs{1}) % this should change to the derivatives folder. Andy wanted to correct this, not sure this works.

if exist('run','var') && run % check if lead dbs should run, not sure this is working correctly

options.dolc = 0;
options.ecog.extractsurface.do = 1;
options.ecog.extractsurface.method = 1;
options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.normalize.do = true;
options.normalize.settings = [];
options.normalize.method = 'ANTs (Avants 2008)';
options.checkreg = false;
options.normalize.refine = 0;
options.coregmr.check = 0;
options.coregmr.method = 'SPM (Friston 2007)';
options.coregmr.do = 1;
options.overwriteapproved = 0;
options.coregct.do = true;
options.coregct.method = 'ANTs (Avants 2008)';
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
options.atl.can = true;
options.atl.pt = 0;
options.atl.ptnative = false;
options.native = 0;
options.d2.col_overlay = 1;
options.d2.con_overlay = 1;
options.d2.con_color = [1 1 1];
options.d2.lab_overlay = 0;
options.d2.bbsize = 50;
options.d2.backdrop = 'MNI152NLin2009bAsym T1 (Fonov 2011)';
options.d2.fid_overlay = 1;
options.d2.write = false;
options.d2.atlasopacity = 0.15;
options.d2.writeatlases = 1;
options.refinelocalization = false;
options.scrf.do = 1;
options.scrf.mask = 'Coarse mask (Schönecker 2008)';
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
options.writeoutpm = 0;
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
