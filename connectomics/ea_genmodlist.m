function [modlist, type] = ea_genmodlist(directory,selectedparc,options)
cnt=1;
modlist=cell(0);
type=[];

% check for canonical fiber sets
fdfibs=dir(ea_getconnectomebase('dmri'));
for fdf=1:length(fdfibs)
    if fdfibs(fdf).isdir && ~strcmp(fdfibs(fdf).name(1),'.')
        [~,fn]=fileparts(fdfibs(fdf).name);
        modlist{cnt}=fn;
        type(cnt)=1;
        cnt=cnt+1;
    end
end

fc=dir(ea_getconnectomebase('fmri'));
for fdf=1:length(fc)
    if fc(fdf).isdir && ~strcmp(fc(fdf).name(1),'.') && ...
       exist([ea_getconnectomebase('fmri'),fc(fdf).name,filesep,'dataset_info.mat'], 'file')
            d=load([ea_getconnectomebase('fmri'),fc(fdf).name,filesep,'dataset_info.mat']);
            [~,fn]=fileparts(fc(fdf).name);
            for ds=1:length(d.dataset.subsets)
                modlist{cnt}=[fn,' > ',d.dataset.subsets(ds).name];
                type(cnt)=2;
                cnt=cnt+1;
            end
    end
end

% patientspecific part
if exist('directory','var') && ~isempty(directory)
    % check if pat-specific fibertracts are present:
    if exist([directory,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized],'file')
        modlist{cnt}='Patient''s fiber tracts';
        type(cnt)=1;
        cnt=cnt+1;
    end

    % fMRI - parcellations:
    % check if rest_tc are present:
    if exist([directory,'connectomics',filesep,selectedparc,filesep,'rest_tc.mat'],'file')
        modlist{cnt}='Patient''s fMRI time courses';
        type(cnt)=2;
        cnt=cnt+1;
    end

    % fMRI - raw files:
    ffis=dir([directory, options.prefs.rest_searchstring]);
    for ff=1:length(ffis)
        [~, restfname] = fileparts(ffis(ff).name);
        modlist{cnt} = ['Patient''s fMRI - ', restfname];
        type(cnt)=2;
        cnt=cnt+1;
    end

end

% check for already processed maps in case of predict module
if exist('options','var')
    if isfield(options,'predict')
        stimdir=[directory,'stimulations',filesep,ea_nt(options),options.predict.stimulation,filesep];
        dfo=dir(stimdir);
        for fo=1:length(dfo)
            if dfo(fo).isdir && ~strcmp(dfo(fo).name(1),'.')
                if ~isempty(dir([stimdir,dfo(fo).name,filesep,'vat_seed_*func_seed_AvgR.nii']))
                    type(cnt)=2; % fMRI result
                    modlist{cnt}=['Precomputed: ',dfo(fo).name];
                    cnt=cnt+1;
                elseif ~isempty(dir([stimdir,dfo(fo).name,filesep,'vat_seed_*struc_seed.nii']))
                    type(cnt)=1; % dMRI result
                    modlist{cnt}=['Precomputed: ',dfo(fo).name];
                    cnt=cnt+1;
                end
            end
        end
    end
end
