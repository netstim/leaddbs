function [modlist, type] = ea_genmodlist(directory,selectedparc,options,conntype)
cnt=1;
modlist=cell(0);
type=[];

% Check 'fmri', 'dmri' or 'both' connectomes
if ~exist('conntype', 'var')
    conntype = 'both';
end

switch lower(conntype)
    case 'both'
        checkdmri = 1;
        checkfmri = 1;
        checkdmri_mt = 1;
    case 'dmri'
        checkdmri = 1;
        checkfmri = 0;
        checkdmri_mt = 0;
    case 'fmri'
        checkdmri = 0;
        checkfmri = 1;
        checkdmri_mt = 0;
    case 'dmri_multitract'
        checkdmri = 0;
        checkdmri_mt = 1;
        checkfmri = 0;

end

% check for canonical fiber sets
if checkdmri
    fdfibs=dir(ea_getconnectomebase('dmri'));
    for fdf=1:length(fdfibs)
        if fdfibs(fdf).isdir && ~strcmp(fdfibs(fdf).name(1),'.')
            [~,fn]=fileparts(fdfibs(fdf).name);
            modlist{cnt}=fn;
            type(cnt)=1;
            cnt=cnt+1;
        end
    end
end

if checkfmri
    fc=dir(ea_getconnectomebase('fmri'));
    for fdf=1:length(fc)
        if fc(fdf).isdir && ~strcmp(fc(fdf).name(1),'.') && ...
           isfile([ea_getconnectomebase('fmri'), fc(fdf).name, filesep, 'dataset_info.json'])
                dataset=loadjson([ea_getconnectomebase('fmri'), fc(fdf).name, filesep, 'dataset_info.json']);
                [~,fn]=fileparts(fc(fdf).name);
                for ds=1:length(dataset.subsets)
                    modlist{cnt}=[fn,' > ',dataset.subsets{ds}.name];
                    type(cnt)=2;
                    cnt=cnt+1;
                end
        end
    end
end
if checkdmri_mt
    fdfibs=dir(ea_getconnectomebase('dMRI_MultiTract'));
    for fdf=1:length(fdfibs)
        if fdfibs(fdf).isdir && ~strcmp(fdfibs(fdf).name(1),'.')
            [~,fn]=fileparts(fdfibs(fdf).name);
            modlist{cnt}=fn;
            type(cnt)=3;
            cnt=cnt+1;
        end
    end
end
% patientspecific part
if exist('directory','var') && ~isempty(directory)
    % check if pat-specific fibertracts are present:
    if checkdmri
        if isfile(fullfile(directory, 'connectomes', 'dMRI', options.prefs.FTR_normalized))
            modlist{cnt}='Patient''s fiber tracts';
            type(cnt)=1;
            cnt=cnt+1;
        end
    end

    % fMRI - parcellations:
    if checkfmri
        % check if rest_tc are present:
        if isfile(fullfile(directory, 'connectomics', selectedparc, 'rest_tc.mat'))
            modlist{cnt}='Patient''s fMRI time courses';
            type(cnt)=2;
            cnt=cnt+1;
        end

        % fMRI - raw files:
        ffis=dir(fullfile(directory, options.prefs.rest_searchstring));
        for ff=1:length(ffis)
            [~, restfname] = fileparts(ffis(ff).name);
            modlist{cnt} = ['Patient''s fMRI - ', restfname];
            type(cnt)=2;
            cnt=cnt+1;
        end
    end
    
    if checkdmri_mt
        if isfile(fullfile(directory, 'connectomes', 'dMRI_MultiTract', options.prefs.FTR_normalized))
            modlist{cnt}='Patient''s fiber tracts';
            type(cnt)=3;
            cnt=cnt+1;
        end
    end
end

% check for already processed maps in case of predict module
if exist('options','var')
    if isfield(options,'predict')
        stimdir = fullfile(directory, 'stimulations', ea_nt(options), options.predict.stimulation);
        dfo=dir(stimdir);
        for fo=1:length(dfo)
            if dfo(fo).isdir && ~strcmp(dfo(fo).name(1),'.')
                if checkfmri && ~isempty(dir(fullfile(stimdir, dfo(fo).name, 'vat_seed_*func_seed_AvgR.nii')))
                    type(cnt)=2; % fMRI result
                    modlist{cnt}=['Precomputed: ',dfo(fo).name];
                    cnt=cnt+1;
                elseif checkdmri && ~isempty(dir(fullfile(stimdir, dfo(fo).name, 'vat_seed_*struc_seed.nii')))
                    type(cnt)=1; % dMRI result
                    modlist{cnt}=['Precomputed: ',dfo(fo).name];
                    cnt=cnt+1;
                end
            end
        end
    end
end
