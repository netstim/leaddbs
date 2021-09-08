function ea_compute_scrf(handles)

options=getappdata(handles.scrf,'options');
directory=getappdata(handles.scrf,'directory');

options.coregmr.method='ANTs';

if ~handles.mask0.Value
    if ~exist([directory,'scrf',filesep,'bgmsk.nii'],'file')
        ea_addtsmask(options,1);
        for msk=1:2
            if msk==2
                btts='2';
            else
                btts='';
            end
            if ~exist([directory,'scrf',filesep],'dir')
                mkdir([directory,'scrf',filesep]);
            end
            movefile([directory,'bgmsk',btts,'.nii'],[directory,'scrf',filesep,'bgmsk',btts,'.nii']);
            ea_conformspaceto([directory,'scrf',filesep,options.prefs.prenii_unnormalized],[directory,'scrf',filesep,'bgmsk',btts,'.nii'],0);
        end
    end
    for msk=1:2
        if msk==2
            btts='2';
        else
            btts='';
        end
        msks{msk}=[directory,'scrf',filesep,'bgmsk',btts,'.nii'];
    end
end
otherfiles={};

if handles.mask1.Value
   msks=msks(1); % only use first mask.
end
if ~exist('msks','var')
    msks={};
end

ea_backuprestore([directory,'scrf',filesep,'movim.nii']);
ea_coregimages(options,[directory,'scrf',filesep,'movim.nii'],[directory,'scrf',filesep,options.prefs.prenii_unnormalized],...
    [directory,'scrf',filesep,'scrfmovim.nii'],otherfiles,1,msks);

%movefile([directory,'scrf',filesep,'movim.nii'],[directory,'scrf',filesep,'scrfmovim.nii']);
movefile([directory,'scrf',filesep,'raw_movim.nii'],[directory,'scrf',filesep,'movim.nii']);

movefile([directory,'scrf',filesep,'movim2',ea_stripext(options.prefs.prenii_unnormalized),'_ants1.mat'],[directory,'scrf',filesep,'scrf_instore.mat']);
delete([directory,'scrf',filesep,ea_stripext(options.prefs.prenii_unnormalized),'2movim','_ants1.mat']);
ea_refreshscrf(options,handles,directory);
