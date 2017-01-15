function ea_applynormtofile_menu(~,~,handles,useinverse)

[options.root,options.patientname]=fileparts(get(handles.patdir_choosebox,'String'));
options.root=[options.root,filesep];
options.earoot=ea_getearoot;
options.prefs=ea_prefs(options.patientname);

[fis,path]=uigetfile({'*.nii';'*.nii.gz'},'Choose files to apply deformation to...',[options.root,options.patientname],'Multiselect','on');

if ~ischar(fis) && ~iscell(fis)
    if ~fis
        return
    end
end
if ischar(fis)
    fis={fis};
end

for fi=1:length(fis) 
   to{fi}=[path,'w',fis{fi}];
   from{fi}=[path,fis{fi}];
end

ea_apply_normalization_tofile(options,from,to,[options.root,options.patientname,filesep],useinverse);

