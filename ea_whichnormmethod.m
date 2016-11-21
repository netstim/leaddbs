function [whichnormmethod,template]=ea_whichnormmethod(directory)
try
    load(fullfile(directory,'ea_normmethod_applied'));
    cnt=0;
    while 1
        whichnormmethod=norm_method_applied{end-cnt};
        switch whichnormmethod
            case {'ea_normalize_apply_normalization','ea_normalize_reslicepretra'}
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
        template=[leaddir,'templates',filesep,'dartel',filesep,'shootmni_6.nii'];
    case 'ea_normalize_spmdartel'
        template=[leaddir,'templates',filesep,'dartel',filesep,'dartelmni_6.nii'];
    case 'ea_normalize_spmnewseg'
        template=[leaddir,'templates',filesep,'TPM_2009b.nii'];
    otherwise
        template=[leaddir,'templates',filesep,'mni_hires.nii'];
end
