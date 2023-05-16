function cs_fmri_conseed_seed_tc(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask)

tic

% if ~isdeployed
%     addpath(genpath('/autofs/cluster/nimlab/connectomes/software/lead_dbs'));
%     addpath('/autofs/cluster/nimlab/connectomes/software/spm12');
% end

if ~exist('writeoutsinglefiles','var')
    writeoutsinglefiles=0;
else
    if ischar(writeoutsinglefiles)
        writeoutsinglefiles=str2double(writeoutsinglefiles);
    end
end


if ~exist('dfold','var')
    dfold=''; % assume all data needed is stored here.
else
    if ~strcmp(dfold(end),filesep)
        dfold=[dfold,filesep];
    end
end

disp(['Connectome dataset: ',cname,'.']);
if contains(cname, '>')
    subset = regexprep(cname, '.*> *', '');
    cname = regexprep(cname, ' *>.*', '');
end

if exist('subset', 'var')
    connLabel = ea_getConnLabel(cname, subset);
else
    connLabel = ea_getConnLabel(cname, 'FullSet');
end

prefs=ea_prefs;
dfoldsurf=[dfold,'fMRI',filesep,cname,filesep,'surf',filesep];
dfoldvol=[dfold,'fMRI',filesep,cname,filesep,'vol',filesep]; % expand to /vol subdir.

dataset = load([dfold,'fMRI',filesep,cname,filesep,'dataset_volsurf.mat']);
dataset.connLabel = connLabel;
dinfo = loadjson([dfold,'fMRI',filesep,cname,filesep,'dataset_info.json']);
dataset.type = dinfo.type;
dataset.subsets = dinfo.subsets;

if exist('outputmask','var')
    if ~isempty(outputmask)
        omask=ea_load_nii(outputmask);
        omaskidx=find(omask.img(:));
        [~,maskuseidx]=ismember(omaskidx,dataset.vol.outidx);
    else
        omaskidx=dataset.vol.outidx;
        maskuseidx=1:length(dataset.vol.outidx);
    end
else
    omaskidx=dataset.vol.outidx; % use all.
    maskuseidx=1:length(dataset.vol.outidx);
end

if ~exist('outputfolder', 'var')
    outputfolder = '';
elseif ~isempty(outputfolder) && ~isfolder(outputfolder)
    mkdir(outputfolder);
end

if strcmp(sfile{1}(end-2:end),'.gz')
    %gunzip(sfile)
    %sfile=sfile(1:end-3);
    usegzip=1;
else
    usegzip=0;
end

for s=1:size(sfile,1)
    if size(sfile(s,:),2)>1
        dealingwithsurface=1;
    else
        dealingwithsurface=0;
    end
    for lr=1:size(sfile(s,:),2)
        if exist(ea_niigz(sfile{s,lr}),'file')
            seed{s,lr}=ea_load_nii(ea_niigz(sfile{s,lr}));
        else
            if size(sfile(s,:),2)==1
                ea_error(['File ',ea_niigz(sfile{s,lr}),' does not exist.']);
            end
            switch lr
                case 1
                    sidec='l';
                case 2
                    sidec='r';
            end
            seed{s,lr}=dataset.surf.(sidec).space; % supply with empty space
            seed{s,lr}.fname='';
            seed{s,lr}.img(:)=0;
        end
        if ~isequal(seed{s,lr}.mat,dataset.vol.space.mat) && (~dealingwithsurface)
            oseedfname=seed{s,lr}.fname;

            try
                seed{s,lr}=ea_conformseedtofmri(dataset,seed{s,lr});
            catch
                keyboard
            end
            seed{s,lr}.fname=oseedfname; % restore original filename if even unneccessary at present.
        end

        [~, seedfn{s,lr}]=fileparts(sfile{s,lr});
        if dealingwithsurface
            sweights=seed{s,lr}.img(:);
        else
            sweights=seed{s,lr}.img(dataset.vol.outidx);
        end
        sweights(isnan(sweights))=0;
        sweights(isinf(sweights))=0; %

        sweights(abs(sweights)<0.0001)=0;
        sweights=double(sweights);

        try
            options=evalin('caller','options');
        end
        if exist('options','var')
            if strcmp(options.lcm.seeddef,'parcellation')
                sweights=round(sweights);
            end
        end
        % assure sum of sweights is 1
        %sweights(logical(sweights))=sweights(logical(sweights))/abs(sum(sweights(logical(sweights))));
        sweightmx=repmat(sweights,1,1);

        sweightidx{s,lr}=find(sweights);
        sweightidxmx{s,lr}=double(sweightmx(sweightidx{s,lr},:));
    end
end

numseed=s;
try
    options=evalin('caller','options');
end
if exist('options','var')
    if strcmp(options.lcm.seeddef,'parcellation') % expand seeds to define
        ea_error('Command not supported for parcellation as input.');
    end
end

disp([num2str(numseed),' seeds, command = ',cmd,'.']);

pixdim=length(dataset.vol.outidx);

numsub=length(dataset.vol.subIDs);

if ~exist('subset','var') % use all subjects
    usesubjects=1:numsub;
else
    for ds=1:length(dataset.subsets)
        if strcmp(subset,dataset.subsets{ds}.name)
            usesubjects=dataset.subsets{ds}.subs;
            break
        end
    end
    numsub=length(usesubjects);
end

% init vars:
for s=1:numseed
    fX{s}=nan(length(omaskidx),numsub);
    rhfX{s}=nan(10242,numsub);
    lhfX{s}=nan(10242,numsub);
end
ea_dispercent(0,'Iterating through subjects');
for sub=1:numsub % iterate across subjects
    howmanyruns=ea_cs_dethowmanyruns(dataset,usesubjects(sub));
    r=cell(howmanyruns,1);
    for run=1:howmanyruns % load data
        r{run}=load([dfoldvol,dataset.vol.subIDs{usesubjects(sub)}{run+1}],'gmtc');
        r{run}.gmtc=single(r{run}.gmtc);
        if isfield(dataset,'surf') && prefs.lcm.includesurf
            if ~exist('ls','var')
                % include surface:
                r{run}.ls=load([dfoldsurf,dataset.surf.l.subIDs{usesubjects(sub)}{run+1}]);
                r{run}.rs=load([dfoldsurf,dataset.surf.r.subIDs{usesubjects(sub)}{run+1}]);
                r{run}.ls.gmtc=single(r{run}.ls.gmtc);
                r{run}.rs.gmtc=single(r{run}.rs.gmtc);
            end
        end
    end

    for s=1:numseed
        thiscorr=zeros(length(omaskidx),howmanyruns);
        if isfield(dataset,'surf') && prefs.lcm.includesurf
            lsthiscorr=zeros(10242,howmanyruns);
            rsthiscorr=zeros(10242,howmanyruns);
        end
        for run=1:howmanyruns
            if size(sfile(s,:),2)>1 % dealing with surface seed
                stc=mean([r{run}.ls.gmtc(sweightidx{s,1},:).*repmat(sweightidxmx{s,1},1,size(r{run}.ls.gmtc,2));...
                    r{run}.rs.gmtc(sweightidx{s,2},:).*repmat(sweightidxmx{s,2},1,size(r{run}.ls.gmtc,2))],1); % seed time course
            else % volume seed
                stc=mean(r{run}.gmtc(sweightidx{s},:).*repmat(sweightidxmx{s},1,size(r{run}.gmtc,2)),1); % seed time course
            end
            thiscorr(:,run)=corr(stc',r{run}.gmtc(maskuseidx,:)','type','Pearson');
            if isfield(dataset,'surf') && prefs.lcm.includesurf
                % include surface:
                lsthiscorr(:,run)=corr(stc',r{run}.ls.gmtc','type','Pearson');
                rsthiscorr(:,run)=corr(stc',r{run}.rs.gmtc','type','Pearson');
            end
        end
        fX{s}(:,sub)=mean(thiscorr,2);
        if isfield(dataset,'surf') && prefs.lcm.includesurf
            lhfX{s}(:,sub)=mean(lsthiscorr,2);
            rhfX{s}(:,sub)=mean(rsthiscorr,2);
        end
        if writeoutsinglefiles
            if isfield(dataset,'surf') && prefs.lcm.includesurf
                ea_writeoutsinglefiles(dataset,outputfolder,sfile,s,usesubjects(sub),thiscorr,omaskidx,lsthiscorr,rsthiscorr)
            else
                ea_writeoutsinglefiles(dataset,outputfolder,sfile,s,usesubjects(sub),thiscorr,omaskidx)
            end
        end
    end

    ea_dispercent(sub/numsub);
end
ea_dispercent(1,'end');

stack = dbstack;
if ismember('ea_networkmapping_calcvals', {stack.name})
    isNetworkMappingRun = 1;
else
    isNetworkMappingRun = 0;
end

for s=1:size(seedfn,1) % subtract 1 in case of pmap command
    % export mean
    M=ea_nanmean(fX{s}',1);
    mmap=dataset.vol.space;
    mmap.dt(1) = 16;
    mmap.img(:)=0;
    mmap.img=single(mmap.img);
    mmap.img(omaskidx)=M;

    if ~isempty(outputfolder)
        mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-AvgR_funcmap.nii']);
    else
        if ~isBIDSFileName(sfile{s})
            mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-AvgR_funcmap.nii']);
        else
            mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'AvgR', 'suffix', 'funcmap');
        end
    end

    if ~isNetworkMappingRun
        ea_write_nii(mmap);
        if usegzip
            gzip(mmap.fname);
            delete(mmap.fname);
        end
    end

    % export variance
    M=ea_nanvar(fX{s}');
    mmap=dataset.vol.space;
    mmap.dt(1) = 16;
    mmap.img(:)=0;
    mmap.img=single(mmap.img);
    mmap.img(omaskidx)=M;

    if ~isempty(outputfolder)
        mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-VarR_funcmap.nii']);
    else
        if ~isBIDSFileName(sfile{s})
            mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-VarR_funcmap.nii']);
        else
            mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'VarR', 'suffix', 'funcmap');
        end
    end

    if ~isNetworkMappingRun
        ea_write_nii(mmap);
        if usegzip
            gzip(mmap.fname);
            delete(mmap.fname);
        end
    end

    if isfield(dataset,'surf') && prefs.lcm.includesurf && ~isNetworkMappingRun
        % lh surf
        lM=ea_nanmean(lhfX{s}');
        lmmap=dataset.surf.l.space;
        lmmap.dt(1) = 16;
        lmmap.img=zeros([size(lmmap.img,1),size(lmmap.img,2),size(lmmap.img,3)]);
        lmmap.img=single(lmmap.img);
        lmmap.img(:)=lM(:);

        if ~isempty(outputfolder)
            lmmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_hemi-L_desc-AvgR_funcmapsurf.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                lmmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_hemi-L_desc-AvgR_funcmapsurf.nii']);
            else
                lmmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'hemi', 'L', 'desc', 'AvgR', 'suffix', 'funcmapsurf');
            end
        end

        ea_write_nii(lmmap);
        if usegzip
            gzip(lmmap.fname);
            delete(lmmap.fname);
        end

        % rh surf
        rM=ea_nanmean(rhfX{s}');
        rmmap=dataset.surf.r.space;
        rmmap.dt(1) = 16;
        rmmap.img=zeros([size(rmmap.img,1),size(rmmap.img,2),size(rmmap.img,3)]);
        rmmap.img=single(rmmap.img);
        rmmap.img(:)=rM(:);

        if ~isempty(outputfolder)
            rmmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_hemi-R_desc-AvgR_funcmapsurf.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                rmmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_hemi-R_desc-AvgR_funcmapsurf.nii']);
            else
                rmmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'hemi', 'R', 'desc', 'AvgR', 'suffix', 'funcmapsurf');
            end
        end

        ea_write_nii(rmmap);
        if usegzip
            gzip(rmmap.fname);
            delete(rmmap.fname);
        end
    end

    % fisher-transform:
    fX{s}=atanh(fX{s});
    if isfield(dataset,'surf') && prefs.lcm.includesurf
        lhfX{s}=atanh(lhfX{s});
        rhfX{s}=atanh(rhfX{s});
    end

    % export fz-mean
    M=nanmean(fX{s}');
    mmap=dataset.vol.space;
    mmap.dt(1) = 16;
    mmap.img(:)=0;
    mmap.img=single(mmap.img);
    mmap.img(omaskidx)=M;

    if ~isempty(outputfolder)
        mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-AvgRFz_funcmap.nii']);
    else
        if ~isBIDSFileName(sfile{s})
            mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-AvgRFz_funcmap.nii']);
        else
            mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'AvgRFz', 'suffix', 'funcmap');
        end
    end

    spm_write_vol(mmap,mmap.img);

    if usegzip
        gzip(mmap.fname);
        delete(mmap.fname);
    end

    if isfield(dataset,'surf') && prefs.lcm.includesurf && ~isNetworkMappingRun
        % lh surf
        lM=nanmean(lhfX{s}');
        lmmap=dataset.surf.l.space;
        lmmap.dt(1) = 16;
        lmmap.img=zeros([size(lmmap.img,1),size(lmmap.img,2),size(lmmap.img,3)]);
        lmmap.img=single(lmmap.img);
        lmmap.img(:)=lM(:);

        if ~isempty(outputfolder)
            lmmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_hemi-L_desc-AvgRFz_funcmapsurf.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                lmmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_hemi-L_desc-AvgRFz_funcmapsurf.nii']);
            else
                lmmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'hemi', 'L', 'desc', 'AvgRFz', 'suffix', 'funcmapsurf');
            end
        end

        ea_write_nii(lmmap);
        if usegzip
            gzip(lmmap.fname);
            delete(lmmap.fname);
        end

        % rh surf
        rM=nanmean(rhfX{s}');
        rmmap=dataset.surf.r.space;
        rmmap.dt(1) = 16;
        rmmap.img=zeros([size(rmmap.img,1),size(rmmap.img,2),size(rmmap.img,3)]);
        rmmap.img=single(rmmap.img);
        rmmap.img(:)=rM(:);

        if ~isempty(outputfolder)
            rmmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_hemi-R_desc-AvgRFz_funcmapsurf.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                rmmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_hemi-R_desc-AvgRFz_funcmapsurf.nii']);
            else
                rmmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'hemi', 'R', 'desc', 'AvgRFz', 'suffix', 'funcmapsurf');
            end
        end

        ea_write_nii(rmmap);
        if usegzip
            gzip(rmmap.fname);
            delete(rmmap.fname);
        end
    end

    % export T
    try
        [~,~,~,tstat]=ttest(fX{s}');
        tmap=dataset.vol.space;
        tmap.img(:)=0;
        tmap.dt(1) = 16;
        tmap.img=single(tmap.img);
        tmap.img(omaskidx)=tstat.tstat;

        if ~isempty(outputfolder)
            tmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-T_funcmap.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                tmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-T_funcmap.nii']);
            else
                tmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'T', 'suffix', 'funcmap');
            end
        end

        if ~isNetworkMappingRun
            spm_write_vol(tmap,tmap.img);
            if usegzip
                gzip(tmap.fname);
                delete(tmap.fname);
            end
        end
    catch
        ea_cprintf('CmdWinWarnings', 'Failed to run connectivity map for seed:\n%s\n', sfile{s});
    end

    if isfield(dataset,'surf') && prefs.lcm.includesurf && ~isNetworkMappingRun
        try
            % lh surf
            [~,~,~,ltstat]=ttest(lhfX{s}');
            lmmap=dataset.surf.l.space;
            lmmap.dt(1) = 16;
            lmmap.img=zeros([size(lmmap.img,1),size(lmmap.img,2),size(lmmap.img,3)]);
            lmmap.img=single(lmmap.img);
            lmmap.img(:)=ltstat.tstat(:);

            if ~isempty(outputfolder)
                lmmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_hemi-L_desc-T_funcmapsurf.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    lmmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_hemi-L_desc-T_funcmapsurf.nii']);
                else
                    lmmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'hemi', 'L', 'desc', 'T', 'suffix', 'funcmapsurf');
                end
            end

            ea_write_nii(lmmap);
            if usegzip
                gzip(lmmap.fname);
                delete(lmmap.fname);
            end
        catch
            ea_cprintf('CmdWinWarnings', 'Failed to run connectivity map (lh surf) for seed:\n%s\n', sfile{s});
        end

        try
            % rh surf
            [~,~,~,rtstat]=ttest(rhfX{s}');
            rmmap=dataset.surf.r.space;
            rmmap.dt(1) = 16;
            rmmap.img=zeros([size(rmmap.img,1),size(rmmap.img,2),size(rmmap.img,3)]);
            rmmap.img=single(rmmap.img);
            rmmap.img(:)=rtstat.tstat(:);

            if ~isempty(outputfolder)
                rmmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_hemi-R_desc-T_funcmapsurf.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    rmmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_hemi-R_desc-T_funcmapsurf.nii']);
                else
                    rmmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'hemi', 'R', 'desc', 'T', 'suffix', 'funcmapsurf');
                end
            end

            ea_write_nii(rmmap);
            if usegzip
                gzip(rmmap.fname);
                delete(rmmap.fname);
            end
        catch
            ea_cprintf('CmdWinWarnings', 'Failed to run connectivity map (rh surf) for seed:\n%s\n', sfile{s});
        end
    end
end

toc


function ea_writeoutsinglefiles(dataset,outputfolder,sfile,s,mcfi,thiscorr,omaskidx,lsthiscorr,rsthiscorr)
ccmap=dataset.vol.space;
ccmap.img=single(ccmap.img);

if ~isempty(outputfolder)
    [~, seedfn] = fileparts(sfile{s});
    ccmap.fname = fullfile(outputfolder, [seedfn, '_conn-', dataset.connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_corr.nii']);
else
    if ~isBIDSFileName(sfile{s})
        ccmap.fname = strrep(sfile{s}, '.nii', ['_conn-', dataset.connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_corr.nii']);
    else
        ccmap.fname = setBIDSEntity(sfile{s}, 'conn', dataset.connLabel, 'id', dataset.vol.subIDs{mcfi}{1}, 'suffix', 'corr');
    end
end

ccmap.img(omaskidx)=mean(thiscorr,2);
ccmap.dt(1) = 16;
spm_write_vol(ccmap,ccmap.img);

% surfs, too:
if isfield(dataset,'surf') && exist('lsthiscorr', 'var')
    ccmap=dataset.surf.l.space;
    ccmap.img=single(ccmap.img);

    if ~isempty(outputfolder)
        [~, seedfn] = fileparts(sfile{s});
        ccmap.fname = fullfile(outputfolder, [seedfn, '_conn-', dataset.connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_hemi-L_corrsurf.nii']);
    else
        if ~isBIDSFileName(sfile{s})
            ccmap.fname = strrep(sfile{s}, '.nii', ['_conn-', dataset.connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_hemi-L_corrsurf.nii']);
        else
            ccmap.fname = setBIDSEntity(sfile{s}, 'conn', dataset.connLabel, 'id', dataset.vol.subIDs{mcfi}{1}, 'hemi', 'L', 'suffix', 'corrsurf');
        end
    end

    ccmap.img(:,:,:,2:end)=[];
    ccmap.img(:)=mean(lsthiscorr,2);
    ccmap.dt(1) = 16;
    spm_write_vol(ccmap,ccmap.img);
end

if isfield(dataset,'surf') && exist('rsthiscorr', 'var')
    ccmap=dataset.surf.r.space;
    ccmap.img=single(ccmap.img);
    ccmap.img(:,:,:,2:end)=[];

    if ~isempty(outputfolder)
        [~, seedfn] = fileparts(sfile{s});
        ccmap.fname = fullfile(outputfolder, [seedfn, '_conn-', dataset.connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_hemi-R_corrsurf.nii']);
    else
        if ~isBIDSFileName(sfile{s})
            ccmap.fname = strrep(sfile{s}, '.nii', ['_conn-', dataset.connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_hemi-R_corrsurf.nii']);
        else
            ccmap.fname = setBIDSEntity(sfile{s}, 'conn', dataset.connLabel, 'id', dataset.vol.subIDs{mcfi}{1}, 'hemi', 'R', 'suffix', 'corrsurf');
        end
    end

    ccmap.img(:)=mean(rsthiscorr,2);
    ccmap.dt(1) = 16;
    spm_write_vol(ccmap,ccmap.img);
end


function howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi)
if strcmp(dataset.type,'fMRI_matrix')
    howmanyruns=1;
else
    howmanyruns=length(dataset.vol.subIDs{mcfi})-1;
end


function X=addone(X)
X=[ones(size(X,1),1),X];


function [mat,loaded]=ea_getmat(mat,loaded,idx,chunk,datadir)

rightmat=(idx-1)/chunk;
rightmat=floor(rightmat);
rightmat=rightmat*chunk;
if rightmat==loaded
    return
end

load([datadir,num2str(rightmat),'.mat']);
loaded=rightmat;
