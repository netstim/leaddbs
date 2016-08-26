function [modlist,sf]=ea_genmodlist(directory,selectedparc,options)
cnt=1;
modlist=cell(0);
% patientspecific part
if exist('directory','var')
    % check if pat-specific fibertracts are present:
    if exist([directory,options.prefs.FTR_normalized],'file');
        modlist{cnt}='Patient-specific fiber tracts';
        sf(cnt)=1;
        cnt=cnt+1;
    end
    
    % fMRI:
    % check if _tc are present:
    if exist([directory,'connectomics',filesep,selectedparc,filesep,'rest_tc.mat'],'file');
        modlist{cnt}='rest_tc';
        sf(cnt)=2;
        cnt=cnt+1;
    end
end


% check for canonical fiber sets
fdfibs=dir([ea_getconnectomebase('dmri'),filesep,'*.mat']);
for fdf=1:length(fdfibs)
    [~,fn]=fileparts(fdfibs(fdf).name);
    modlist{cnt}=fn;
    sf(cnt)=1;
    cnt=cnt+1;
end

fc=dir([ea_getconnectomebase('fmri'),filesep,'']);
for fdf=1:length(fc)
    if fc(fdf).isdir && ~strcmp(fc(fdf).name(1),'.')
    [~,fn]=fileparts(fc(fdf).name);
    modlist{cnt}=fn;
    sf(cnt)=2;
    cnt=cnt+1;
    end
end
