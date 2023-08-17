function ea_run_deface(hobj, evt, handles)

answ=questdlg(['Warning: this method only works for the MNI152NLin2009bAsym Space!' sprintf('\n\n') 'Do you want to run brain masking or defacing to anonymize patient data?'], ...
    'De-identify base files', ...
    'Mask','Deface','Cancel', 'Cancel');

if ~strcmp(answ,'Cancel')
    if strcmp( ea_getspace, 'MNI152NLin2009bAsym')      % check if space is actually MNI ICBM 2009b
        dirs=getappdata(handles.leadfigure, 'uipatdir');                 % get current list of patients loaded in the gui
        options = ea_getptopts(dirs{1});
        [pth,fn]=fileparts(options.bids.datasetDir);
        if ~exist(fullfile(pth,[fn,'_deidentified']),'dir')
            ea_mkdir(fullfile(pth,[fn,'_deidentified']));
        end
        outdir=uigetdir(fullfile(pth,[fn,'_deidentified']),'Choose output folder');
        for pt=1:length(dirs)
            options = ea_getptopts(dirs{pt});
            nm=loadjson(options.subj.norm.log.method);
            whichnormmethod=nm.method;
            newID=inputdlg(['Enter Non-linked ID for ',options.subj.subjId],'De-identify',1,{options.subj.subjId});
            switch whichnormmethod                   % before running, check whether patient has been normalized or not
                case ''
                    fprintf('\nDe-idenfitication not possible for patient %s, because normalization was not performed!\nPlease perform normalization first.\n', options.patientname);
                otherwise
                    fprintf('\nRunning de-idenfitication with the option ''%s''...\n', answ);
                    ea_run_anonymize_pt(dirs{pt}, answ, options, outdir,newID{1});               % run deface function with desired methode (deface or brain extract)
                    fprintf('** Process Done.\n')
            end

        end
    else
        fprintf('Space is %s, please switch to MNI152NLin2009bAsym space to enable de-identification of patient data!\n', ea_getspace)
    end
end

function  ea_run_anonymize_pt(directory, method, options, outdir, newID)

if length(newID)<4
    newID=['sub-',newID];
else
    if ~startsWith(newID,'sub-')
        newID=['sub-',newID];
    end
end

ea_mkdir(fullfile(outdir,'derivatives','leaddbs'));
[pth,oldID]=fileparts(directory);
if exist(fullfile(outdir,'derivatives','leaddbs',newID),'dir')
    ea_error([fullfile(outdir,'derivatives','leaddbs',newID),' already exists. Please clear the output destination.']);
end
ea_mkdir(fullfile(outdir,'derivatives','leaddbs',newID));
% copy files to new folder
newdirectory=fullfile(outdir,'derivatives','leaddbs',[newID]);
list=ea_bidsfiles2list(options);

for file=1:length(list)
    thisinfile=list{file};
    thisoutfile=strrep(thisinfile,oldID,newID);
    thisoutfile=strrep(thisoutfile,options.bids.datasetDir,outdir);
    [pth,fn]=fileparts(thisoutfile);
    ea_mkdir(pth);
    if exist(thisinfile,'file')
        copyfile(thisinfile,thisoutfile);
    end
end

% switch pointers to new folder
 options=ea_getptopts(newdirectory);
directory=newdirectory;

switch method
    % simple brain masking using the MNI ICBM 2009b brain mask
    case 'Mask'
        mask_file = fullfile(options.earoot, 'templates', 'space', ea_getspace, 'brainmask.nii.gz');
        % defacing by defaced MNI ICBM 2009b brain mask
    case 'Deface'
        mask_file = fullfile(options.earoot, 'templates', 'space', ea_getspace, 'brainmask_defaced.nii.gz');
end
mask_file_pt_space = fullfile(ea_getleadtempdir,['patient_mask',ea_generate_uuid,'.nii.gz']);

% normalize MNI brainmask to patient space
ea_apply_normalization_tofile(options, mask_file, mask_file_pt_space, 1, 0, options.subj.coreg.anat.preop.(options.subj.AnchorModality));

% now go through all the files and multiply them with the mask to
% create brain extracted images

% PREOP subject space
fn=fieldnames(options.subj.coreg.anat.preop);
for fi = 1:length(fn)
    thisfile=fullfile(options.subj.coreg.anat.preop.(fn{fi}));
    if exist(thisfile, 'file')
        applyMask(thisfile, thisfile, mask_file_pt_space);
    end
end

% POSTOP subject space
fn=fieldnames(options.subj.coreg.anat.postop);
for fi = 1:length(fn)
    thisfile=options.subj.coreg.anat.postop.(fn{fi});
    if exist(thisfile, 'file')
        applyMask(thisfile, thisfile, mask_file_pt_space);
    end
end

% PREOP template space
fn=fieldnames(options.subj.norm.anat.preop);
for fi = 1:length(fn)
    thisfile=options.subj.norm.anat.preop.(fn{fi});
    if exist(thisfile, 'file')
        applyMask(thisfile, thisfile, mask_file);
    end
end

% POSTOP template space
fn=fieldnames(options.subj.norm.anat.postop);
for fi = 1:length(fn)
    thisfile=options.subj.norm.anat.postop.(fn{fi});
    if exist(thisfile, 'file')
        applyMask(thisfile, thisfile, mask_file);
    end
end

% apply bet to raw files:
fn=fieldnames(options.subj.preopAnat);
for fi = 1:length(fn)
    thisfile=options.subj.preopAnat.(fn{fi});
    ea_bet(thisfile.raw);
    ea_bet(thisfile.preproc);
end
fn=fieldnames(options.subj.postopAnat);
for fi = 1:length(fn)
    thisfile=options.subj.postopAnat.(fn{fi});
    ea_bet(thisfile.raw);
    ea_bet(thisfile.preproc);
end

function applyMask(infilename, outfilename, mask_filename)
mask_pt_space = ea_load_nii(mask_filename);         % load mask .nii
mask_pt_space.img=logical(mask_pt_space.img);    % convert to logical

anatfile = ea_load_nii(infilename);
outfile = anatfile;
outfile.img(~logical(mask_pt_space.img)) =0;
outfile.fname=outfilename;
ea_write_nii(outfile);
fprintf('Written de-identified image: %s\n', outfilename);









