function [whichnormmethod,tempfile]=ea_whichnormmethod(directory)
try
load([directory,'ea_normmethod_applied']);
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

leaddir=[ea_getearoot];

switch whichnormmethod
    case 'ea_normalize_spmdartel' % use dartel MNI template
        tempfile=[leaddir,'templates',filesep,'dartel',filesep,'dartelmni_6.nii'];
    case 'ea_normalize_spmnewseg'
        tempfile=[leaddir,'templates',filesep,'TPM.nii'];
    otherwise % use mni_hires
        tempfile=[leaddir,'templates',filesep,'mni_hires.nii'];
end
