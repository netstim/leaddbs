function ea_run_deface(hobj, evt, handles)

answ=questdlg('Do you want to run brain extraction or defacing?', ...
    'Deface base files', ...
    'Extract','Deface','Cancel', 'Cancel');

if ~strcmp(answ,'Cancel')   
    dirs=getappdata(handles.leadfigure,'uipatdir');     % get current list of patients loaded in the gui
    
    for pt=1:length(dirs)
        ea_run_deface_pt(dirs{pt}, answ);               % run deface function with desired methode (deface or brain extract)
    end
end

function ea_run_deface_pt(directory, method)

options=ea_getptopts(directory);
options.coregmr.method='ANTs'; % hardcode for now.
[~,presentfiles] = ea_assignpretra(options);
presentfiles=stripstn(presentfiles);

presentfiles=ea_addpostops(presentfiles,options);

if ~strcmp(directory(end),filesep)
    directory=[directory,filesep];
end

tempdir=ea_getleadtempdir;
copyfile([ea_space,'mask.nii.gz'],[tempdir,'mask.nii.gz']);
gunzip([tempdir,'mask.nii.gz']);
delete([tempdir,'mask.nii.gz']);
ea_apply_normalization_tofile(options,{[tempdir,'mask.nii']},{[directory,'brainmask.nii']},directory,1,0,[directory,presentfiles{1}]);

mask=ea_load_nii([directory,'brainmask.nii']);
for fi=1:length(presentfiles)
    
    
    % also deface raw pre-/postop files
    if exist([directory,'raw_',presentfiles{fi}],'file');
        defaceraw([directory,presentfiles{fi}],[directory,'raw_',presentfiles{fi}],directory,options);
    end
    
    % also deface raw postop-ct file
    if exist([directory,presentfiles{fi}(2:end)],'file') && strcmp(presentfiles{fi}(2:end),options.prefs.rawctnii_unnormalized);
        defaceraw([directory,presentfiles{fi}],[directory,presentfiles{fi}(2:end)],directory,options);
    end

    
    if exist([directory,presentfiles{fi}],'file');
        thisanat=ea_load_nii([directory,presentfiles{fi}]);
        thisanat.img(~logical(mask.img))=0;
        ea_write_nii(thisanat);
    end
end

function defaceraw(registered,rawfile,directory,options)
 % going back by backregistering here. We should change this to export
    % .mats by default when registering raw to anchor modality in the first
    % place, in the future.

    copyfile([directory,'brainmask.nii'],[directory,'rbrainmask.nii']); % create brainmask for raw file copy.
    ea_coreg2images(options,registered,rawfile,[directory,'tmp.nii'],{[directory,'rbrainmask.nii']},0,[],0);
    
    rmask=ea_load_nii([directory,'rbrainmask.nii']);
    rmask.img=logical(rmask.img);
    
    thisanat=ea_load_nii(rawfile);
    thisanat.img(~logical(rmask.img))=0;
    ea_write_nii(thisanat);



function presentfiles=stripstn(presentfiles)

[mem,ix]=ismember(presentfiles,{'anat_STN.nii','anat_RN.nii','anat_GPi.nii','anat_GPe.nii'});
presentfiles=presentfiles(~mem);

function presentfiles=ea_addpostops(presentfiles,options)

switch options.modality
    case 1 % MR
        presentfiles{end+1}=options.prefs.tranii_unnormalized;
        presentfiles{end+1}=options.prefs.cornii_unnormalized;
        presentfiles{end+1}=options.prefs.sagnii_unnormalized;
    case 2 % CT
        presentfiles{end+1}=options.prefs.ctnii_coregistered;
end
                        
                        
                        
                        
                        
                        
                        
                        
                        