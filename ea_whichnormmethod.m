function [whichnormmethod,tempfile]=ea_whichnormmethod(directory)
load([directory,'ea_normmethod_applied']);
cnt=0;
while 1
    whichnormmethod=norm_method_applied{end-cnt};
    switch whichnormmethod
        case 'ea_normalize_apply_normalization'
            cnt=cnt+1;
        otherwise
            break
    end

end

leaddir=[fileparts(which('lead')),filesep];

switch whichnormmethod
    case 'ea_normalize_spmdartel' % use dartel MNI template
        tempfile=[leaddir,'templates',filesep,'dartel',filesep,'dartelmni_6.nii'];
%     case 'ea_normalize_ants'
%         ea_error('ANTs normalization is not supported as of now.');
    otherwise % use mni_hires
        tempfile=[leaddir,'templates',filesep,'mni_hires.nii'];
end