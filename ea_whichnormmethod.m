function [whichnormmethod,template]=ea_whichnormmethod(directory,whatmethod)

if ~strcmp(directory(end),filesep)
   directory=[directory,filesep];
end
whichnormmethod=''; % default exit empty.
if ~exist('whatmethod','var')
    whatmethod='normmethod';
end

switch whatmethod
    case 'normmethod'
        varname='norm_method_applied';

        try
            load(fullfile(directory,['ea_',whatmethod,'_applied']));
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
    case 'coregctmethod'

        varname='coregct_method_applied';
        try
            load(fullfile(directory,['ea_',whatmethod,'_applied']));
            whichnormmethod=eval([varname]);

            if ~iscell(whichnormmethod)
                whichnormmethod={whichnormmethod};
            end
        end

    case 'coregmrmethod'
        varname='coregmr_method_applied';
        try
            load(fullfile(directory,['ea_',whatmethod,'_applied']));
            whichnormmethod=eval([varname]);
            if ~iscell(whichnormmethod)
                whichnormmethod={whichnormmethod};
            end
        end
end

