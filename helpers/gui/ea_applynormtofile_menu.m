function ea_applynormtofile_menu(~, ~, handles, useinverse, untouchedanchor, asoverlay, expdicom, fname, templateresolution, targetfile)
if ~exist('untouchedanchor','var')
    untouchedanchor = 0;
end

if ~exist('templateresolution','var')
    templateresolution = 0;
end

if ~exist("targetfile",'var')
    targetfile=0;
end

if templateresolution
    res = inputdlg('Specify voxel resolution of template space to warp into.','Template resolution',1,{'0.5'});
    templateresolution = str2double(res);
end


if ~exist('expdicom','var')
    expdicom = 0;
end

if ~exist('asoverlay','var')
    asoverlay = 0;
end

if ~iscell(handles)
    uipatdir = getappdata(handles.leadfigure,'uipatdir');
else
    uipatdir = handles; % direct supply of cell string.
end
if ~exist('fname','var')
    fname=0;
end

if ~fname || isempty(fname)
    if useinverse
        defaultPath = ea_space;
    else
        if length(uipatdir) == 1
            defaultPath = uipatdir{1};
        else
            defaultPath = fileparts(uipatdir{1});
        end
    end
    if ismac % macs file open dlg doesnt seem to support .nii/.nii.gz handling from matlab.
        [fromfis, frompath] = uigetfile({'*'}, 'Choose files to apply deformation to...', defaultPath, 'Multiselect', 'on');
    else
        [fromfis, frompath] = uigetfile({'*.nii' 'NIfTI';'*.nii.gz' 'Compressed NIfTI'}, 'Choose files to apply deformation to...', defaultPath, 'Multiselect', 'on');
    end
    if ~ischar(fromfis) && ~iscell(fromfis)
        if ~fromfis
            return
        end
    end
else
    [frompath, fromfis, ext] = fileparts(fname);
    if ~isempty(frompath)
        frompath = fullfile(frompath, filesep);
    else % local file
        frompath = fullfile('.', filesep);
    end
    fromfis = [fromfis, ext];
end

if ischar(fromfis)
    fromfis = {fromfis};
end
space=''; % default blank
if useinverse % from template space to [untouched] achor space
    for pt=1:length(uipatdir)
        options = ea_getptopts(uipatdir{pt});
        presentfiles = fieldnames(options.subj.preopAnat);
        options.coregmr.method = get_coregmr_method;

        from = cell(length(fromfis), 1);
        to = cell(length(fromfis), 1);

        if untouchedanchor
            spaceTag = [options.subj.subjId, 'Native'];
        else
            spaceTag = [options.subj.subjId, 'anchorNative'];
        end

        for i=1:length(fromfis)
        
                if isBIDSFileName(fromfis{i})
                    to{i} = setBIDSEntity(fullfile(frompath, fromfis{i}), 'space', spaceTag);
                else
                    to{i} = strrep(fromfis{i}, '.nii', ['_space-', spaceTag, '.nii']);
                end
           
                if targetfile % overwrite since will be supplied in better resolution
                    if ischar(targetfile)
                        to{i}=targetfile;
                    elseif iscell(targetfile)
                        to{i}=targetfile{i};
                    else
                        if length(fromfis)>1
                            ea_error('Not supported for multiple images, please supply one by one.');
                        end
                        [spacefile,spacepath]=uigetfile({'*.nii' 'NIfTI';'*.nii.gz' 'Compressed NIfTI'},['Specify target space for image number ',num2str(i)]);
                        space=fullfile(spacepath,spacefile);
                    end
                end
            from{i} = fullfile(frompath, fromfis{i});
        end
        if length(from)==1
            interp='auto';
        else
            interp=1;
        end
        ea_apply_normalization_tofile(options, from, to, useinverse, interp, space);

        if untouchedanchor % map from anchor to untouched anchor
            tmp_file = strrep(options.subj.preproc.anat.preop.(presentfiles{1}),'desc-preproc','desc-tmp');
            ea_coregimages(options,...
                options.subj.coreg.anat.preop.(presentfiles{1}),...
                options.subj.preproc.anat.preop.(presentfiles{1}),...
                tmp_file, to,[],[],interp);
            ea_delete(tmp_file);

            if asoverlay
                untouchedanchorImage=ea_load_nii(options.subj.preproc.anat.preop.(presentfiles{1}));
                overlay=ea_load_nii(to{1});
                fused=untouchedanchorImage;
                fused.img(:)=zscore(fused.img(:));
                fused.img=fused.img+overlay.img;
                fused.img=ea_rescale(fused.img);
                fused.img=fused.img*255;
                fused.dt(1) = 2;
                [natpath,natfn,natext]=fileparts(untouchedanchorImage.fname);
                fused.fname=fullfile(natpath,[natfn,'_overlay',natext]);
                ea_write_nii(fused);
            end

            if expdicom
                natpath=fileparts(untouchedanchorImage.fname);
                [filename,pathname]=uigetfile('*.*','Select sample DICOM',[natpath,filesep,'DICOM']);
                dicom_file=fullfile(pathname,filename);
                merged_file=fused.fname;
                newSeriesNumber=100;
                newSeriesDescription='LeadDBS Plan';
                mkdir(fullfile(natpath,'DICOM','LeadDBSExport'));
                outputDirectory=fullfile(natpath,'DICOM','LeadDBSExport');
                mergedImageVolume=1;
                outputImagePosition=2;

                uw_overlay_convert2dicom(dicom_file, merged_file, newSeriesNumber, newSeriesDescription, outputDirectory, mergedImageVolume, outputImagePosition);
            end
        end
    end
    
else % from [untouched] achor space to template space
    options = ea_getptopts(uipatdir{1});
    presentfiles = fieldnames(options.subj.preopAnat);
    options.coregmr.method = get_coregmr_method;

    to = cell(length(fromfis), 1);
    for i=1:length(fromfis)
        if isBIDSFileName(fromfis{i})
            to{i} = setBIDSEntity(fullfile(frompath, fromfis{i}), 'space', ea_getspace);
        else
            to{i} = strrep(fromfis{i}, '.nii', ['_space-', ea_getspace, '.nii']);
        end
        copyfile(fullfile(frompath, fromfis{i}), to{i});
    end

    if untouchedanchor % map from untouched anchor to anchor first
        tmp_file = strrep(options.subj.preproc.anat.preop.(presentfiles{1}),'desc-preproc','desc-tmp');
        ea_coregimages(options,...
            options.subj.preproc.anat.preop.(presentfiles{1}),...
            options.subj.coreg.anat.preop.(presentfiles{1}),...
            tmp_file, to);
        ea_delete(tmp_file);
    end

    if templateresolution
        ea_mkdir([ea_space,'resliced_templates']);
        trstr=num2str(templateresolution);
        trstr=strrep(trstr,'.','_');
        if ~exist([ea_space,'resliced_templates',filesep,trstr,'.nii.gz'],'file')
            copyfile(ea_niigz([ea_space,options.primarytemplate]),[ea_space,'resliced_templates',filesep,trstr,'.nii']);
            ea_reslice_nii([ea_space,'resliced_templates',filesep,trstr,'.nii'],[ea_space,'resliced_templates',filesep,trstr,'.nii'],repmat(templateresolution,1,3));
            nii=ea_load_nii([ea_space,'resliced_templates',filesep,trstr,'.nii']);
            nii.img(:)=0;
            nii.dt(1) = 4;
            ea_write_nii(nii);
            gzip(nii.fname);
            delete(nii.fname);
        end
        refim=[ea_space,'resliced_templates',filesep,trstr,'.nii.gz'];
    else
        refim=ea_niigz([ea_space,options.primarytemplate]);
    end
    if length(from)==1
        interp='auto';
    else
        interp=1;
    end
    ea_apply_normalization_tofile(options, to, to, useinverse, interp, refim);
end


function coregmr_method = get_coregmr_method(handles)
try
    coregmr_method = handles.coregctmethod.String{handles.coregmrmethod.Value};
catch
    coregmr_method = 'SPM (Friston 2007)';
end
