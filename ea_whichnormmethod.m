function [whichnormmethod,template]=ea_whichnormmethod(directory,whatmethod)

if ~strcmp(directory(end),filesep)
   directory=[directory,filesep];
end

whichnormmethod=''; % default exit empty.
if ~exist('whatmethod','var')
    whatmethod='normmethod';
end

options = ea_getptopts(directory);
    
switch whatmethod
    case 'normmethod'

        try
            json = loadjson(options.subj.norm.log.method);
            whichnormmethod = upper(regexp(json.method, '\S+', 'match', 'once'));
        catch
            whichnormmethod='';
        end

        if contains(lower(json.method), 'shoot')
            template=[ea_space([],'dartel'),'shootmni_6.nii'];
        elseif contains(lower(json.method), 'dartel')
            template=[ea_space([],'dartel'),'shootmni_6.nii'];
        elseif contains(lower(json.method), 'newseg')
            template=[ea_space,'TPM.nii'];
        else
            spacedef=ea_getspacedef;
            template=[ea_space,spacedef.templates{1},'.nii'];
        end
        
    case 'coregctmethod'

        json = loadjson(options.subj.coreg.log.method);
        whichnormmethod = json.method.CT;

    case 'coregmrmethod'
        
        json = loadjson(options.subj.coreg.log.method);
        fnames = fieldnames(json.method);
        mrfname = fnames{find(~strcmp(fnames,'CT'),1)};
        whichnormmethod = json.method.(mrfname);
        
end

