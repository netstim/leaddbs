function ea_applynormto_custom_file(uipatdir,  untouchedanchor, fname, interp)
% interp = interpolation to choose (default = 1) - losely following spm nomenclature (0 = nearest neighbor, 1 = trilinear etc)
% untouchedanchor = 0 -> use anchor space, 1 -> use raw anchor space, 'filename.nii' -> use specific file
% fname = output file name
% templateresolution = which isotropic resolution to choose (default = 0.5mm)


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

if ~exist('interp','var')
    interp=1;
end

if ~exist('fname','var') || isempty(fname)
    
        if length(uipatdir) == 1
            defaultPath = uipatdir{1};
        else
            defaultPath = fileparts(uipatdir{1});
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
    fis={fis};
    path={path};
else
    if ~iscell(fname)
        fname={fname};
    end
    [path,fis,ext]=fileparts(fname);
    for f=1:length(path)
        path{f}=[path{f}, filesep];
        fis{f}=[fis{f},ext{f}];
    end
end

[options.root, options.patientname] = fileparts(uipatdir);
options.root = [options.root, filesep];
options.earoot = ea_getearoot;
options.prefs = ea_prefs(options.patientname);
[options,presentfiles] = ea_assignpretra(options);

options.coregmr.method='SPM';


if ischar(fis)
    fis = {fis};
end

to = cell(1, length(fis));
from = cell(1, length(fis));

for fi=1:length(fis)
    to{fi} = [path{fi}, 'gl', fis{fi}];
    from{fi} = [path{fi}, fis{fi}];
end

if untouchedanchor % map from untouched anchor to anchor first
    if ischar(untouchedanchor)
        input=fullfile(options.root, options.patientname, untouchedanchor);
    else
        input=[options.root, options.patientname, filesep, presentfiles{1}];
    end
    uid=ea_generate_uuid;
    copyfile(input,[ea_getleadtempdir,uid,'.nii'])

    for f=1:length(from)
        [pp,ff,ee]=fileparts(from{f});
        copyfrom{f}=fullfile(pp,['r',ff,ee]);
        copyfile(from{f},copyfrom{f});
    end

    ea_coreg2images(options,[ea_getleadtempdir,uid,'.nii'],...
        [options.root, options.patientname, filesep, presentfiles{1}],...
        [ea_getleadtempdir,'r',uid,'.nii'],...
        copyfrom);
    ea_delete([ea_getleadtempdir,uid,'.nii']);
    ea_delete([ea_getleadtempdir,'r',uid,'.nii']);
    from=copyfrom;
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

ea_apply_normalization_tofile(options, from, to, [options.root, options.patientname, filesep], 0, interp, refim);
