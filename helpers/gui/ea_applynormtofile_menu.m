function ea_applynormtofile_menu(~, ~, handles, useinverse, untouchedanchor, asoverlay, expdicom, fname, templateresolution)
if ~exist('untouchedanchor','var')
    untouchedanchor=0;
end
if ~exist('templateresolution','var')
    templateresolution=0;
end
if templateresolution
   res=inputdlg('Specify voxel resolution of template space to warp into.','Template resolution',1,{'0.5'});
   templateresolution=str2double(res);
end

if untouchedanchor
    interp=0;
else
    interp=4;
end

if ~exist('expdicom','var')
    expdicom=0;
end

if ~exist('asoverlay','var')
    asoverlay=0;
end

if ~iscell(handles)
    uipatdir=getappdata(handles.leadfigure,'uipatdir');
else
    uipatdir=handles; % direct supply of cell string.
end

if ~exist('fname','var') || isempty(fname)
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
        [fis, path] = uigetfile({'*'}, 'Choose files to apply deformation to...', defaultPath, 'Multiselect', 'on');
    else
        [fis, path] = uigetfile({'*.nii' 'NIfTI';'*.nii.gz' 'Compressed NIfTI'}, 'Choose files to apply deformation to...', defaultPath, 'Multiselect', 'on');
    end
    if ~ischar(fis) && ~iscell(fis)
        if ~fis
            return
        end
    end
else
    [path,fis,ext]=fileparts(fname);
    if ~isempty(path)
        path=[path, filesep];
    else % local file
        path=['.', filesep];
    end
    fis=[fis,ext];
end

if useinverse % from template space to [untouched] achor space
    for pt=1:length(uipatdir)
        [options.root, options.patientname] = fileparts(uipatdir{pt});
        options.root = [options.root, filesep];
        options.earoot = ea_getearoot;
        options.prefs = ea_prefs(options.patientname);
        [options,presentfiles] = ea_assignpretra(options);
        try
            options.coregmr.method=get(handles.coregmrmethod,'String');
            options.coregmr.method=options.coregmr.method{get(handles.coregmrmethod,'Value')};
        catch
            options.coregmr.method='SPM (Friston 2007)';
        end

        if ischar(fis)
            fis = {fis};
        end

        to = cell(1, length(fis));
        from = cell(1, length(fis));

        for fi=1:length(fis)
            to{fi} = [uipatdir{pt}, filesep, 'w', fis{fi}];
            from{fi} = [path, fis{fi}];
        end

        ea_apply_normalization_tofile(options, from, to, [options.root, options.patientname, filesep], useinverse, interp);

        if untouchedanchor % map from anchor to untouched anchor
            ea_coregimages(options,[options.root, options.patientname, filesep, presentfiles{1}],...
                [options.root, options.patientname, filesep, 'raw_', presentfiles{1}],...
                [options.root, options.patientname, filesep, 'rraw_', presentfiles{1}],...
                to,[],[],interp);
            ea_delete([options.root, options.patientname, filesep, 'rraw_', presentfiles{1}]);

            if asoverlay
                untouchedanchorImage=ea_load_nii([options.root, options.patientname, filesep, 'raw_', presentfiles{1}]);
                overlay=ea_load_nii(to{1});
                fused=untouchedanchorImage;
                fused.img(:)=zscore(fused.img(:));
                fused.img=fused.img+overlay.img;
                fused.img=ea_rescale(fused.img);
                fused.img=fused.img*255;
                fused.dt=[2,0];
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
    [options.root, options.patientname] = fileparts(fileparts(path));
    options.root = [options.root, filesep];
    options.earoot = ea_getearoot;
    options.prefs = ea_prefs(options.patientname);
    [options,presentfiles] = ea_assignpretra(options);
    try
        options.coregmr.method=get(handles.coregmrmethod,'String');
        options.coregmr.method=options.coregmr.method{get(handles.coregmrmethod,'Value')};
    catch
        options.coregmr.method='SPM (Friston 2007)';
    end

    if ischar(fis)
        fis = {fis};
    end

    to = cell(1, length(fis));
    from = cell(1, length(fis));

    for fi=1:length(fis)
        to{fi} = [path, 'gl', fis{fi}];
        from{fi} = [path, fis{fi}];
    end

    if untouchedanchor % map from untouched anchor to anchor first
        ea_coregimages(options,[options.root, options.patientname, filesep, 'raw_', presentfiles{1}],...
            [options.root, options.patientname, filesep, presentfiles{1}],...
            [options.root, options.patientname, filesep, 'uraw_', presentfiles{1}],...
            from);
        ea_delete([options.root, options.patientname, filesep, 'uraw_', presentfiles{1}]);
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
            nii.dt=[4,0];
            ea_write_nii(nii);
            gzip(nii.fname);
            delete(nii.fname);
        end
        refim=[ea_space,'resliced_templates',filesep,trstr,'.nii.gz'];
    else
        refim=ea_niigz([ea_space,options.primarytemplate]);
    end

    ea_apply_normalization_tofile(options, from, to, [options.root, options.patientname, filesep], useinverse, interp, refim);
end
