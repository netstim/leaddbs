function [whichnormmethod,template]=ea_coregmethod(directory)
try
    load(fullfile(directory,'ea_normmethod_applied'));
    cnt=0;
    while 1
        whichnormmethod=norm_method_applied{end-cnt};
        switch whichnormmethod
            case {'ea_normalize_apply_normalization'}
                cnt=cnt+1;
            otherwise
                break
        end
    end
catch
    whichnormmethod='';
end

leaddir=ea_getearoot;
switch whichnormmethod
    case 'ea_normalize_spmshoot'
        template=[ea_space([],'dartel'),'shootmni_6.nii'];
    case 'ea_normalize_spmdartel'
        template=[ea_space([],'dartel'),'dartelmni_6.nii'];
    case 'ea_normalize_spmnewseg'
        template=[ea_space,'TPM.nii'];
    otherwise
        options.prefs=ea_prefs('');
        spacedef=ea_getspacedef;
        template=[ea_space,spacedef.templates{1},'.nii'];
end
