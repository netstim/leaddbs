function ea_applynormtofile_menu(~, ~, handles, useinverse, untouchedanchor, asoverlay,expdicom)
if ~exist('untouchedanchor','var')
    untouchedanchor=0;
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

[options.root, options.patientname] = fileparts(get(handles.patdir_choosebox, 'String'));
options.root = [options.root, filesep];
options.earoot = ea_getearoot;
options.prefs = ea_prefs(options.patientname);
[options,presentfiles] = ea_assignpretra(options);
options.coregmr.method='SPM';
[fis, path] = uigetfile({'*.nii';'*.nii.gz'}, 'Choose files to apply deformation to...', [options.root, options.patientname], 'Multiselect', 'on');

if ~ischar(fis) && ~iscell(fis)
    if ~fis
        return
    end
end
if ischar(fis)
    fis = {fis};
end

to = cell(1, length(fis));
from = cell(1, length(fis));

if untouchedanchor && ~useinverse % need to map from untouched anchor anchor first
    keyboard
    ea_coreg2images_generic(options,[options.root, options.patientname, filesep, 'raw_',presentfiles{1}],...
        [options.root, options.patientname, filesep, presentfiles{1}],...
        [options.root, options.patientname, filesep, 'uraw_',presentfiles{1}],...
        from);
end

for fi=1:length(fis)
    to{fi} = [path, 'gl', fis{fi}];
    from{fi} = [path, fis{fi}];
end

ea_apply_normalization_tofile(options, from, to, [options.root, options.patientname, filesep], useinverse, interp);

if untouchedanchor && useinverse % need to map from anchor to untouched (raw) anchor
    ea_coreg2images_generic(options,[options.root, options.patientname, filesep, presentfiles{1}],...
        [options.root, options.patientname, filesep, 'raw_',presentfiles{1}],...
        [options.root, options.patientname, filesep, 'rraw_',presentfiles{1}],...
        to);
    if asoverlay
        keyboard
        
    end
    
    if expdicom
        uigetfile(
        dicom_file='/Users/leaddbs/Downloads/uw_overlay/t1_mprage_tra_p2_iso_8/IM-0004-0192.dcm';
        merged_file='temp.img';
        newSeriesNumber=100;
        newSeriesDescription='Lead-DBS Plan';
        outputDirectory='/Users/leaddbs/Downloads/uw_overlay/';
        mergedImageVolume=1;
        outputImagePosition=2;
        
        uw_overlay_convert2dicom(dicom_file, merged_file, newSeriesNumber, newSeriesDescription, outputDirectory, mergedImageVolume, outputImagePosition);
 
        
    end
end

