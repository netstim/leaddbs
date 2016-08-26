function modlist=ea_genmodlist(directory,selectedparc,options)
cnt=1;
modlist=cell(0);
% check if pat-specific fibertracts are present:
if exist([directory,options.prefs.FTR_normalized],'file');
    modlist{cnt}='Patient-specific fiber tracts';
    cnt=cnt+1;
end
% check for canonical fiber sets
fdfibs=dir([options.earoot,ea_getconnectomebase('dmri'),filesep,'*.mat']);
for fdf=1:length(fdfibs)
    [~,fn]=fileparts(fdfibs(fdf).name);
    modlist{cnt}=fn;
    cnt=cnt+1;
end

% fMRI:
% check if _tc are present:
if exist([directory,'connectomics',filesep,selectedparc,filesep,'rest_tc.mat'],'file');
    modlist{cnt}='rest_tc';
end