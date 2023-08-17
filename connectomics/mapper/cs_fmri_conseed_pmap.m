function cs_fmri_conseed_pmap(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask)

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

dfoldvol=[dfold,'fMRI',filesep,cname,filesep,'vol',filesep]; % expand to vol subdir.

dataset = load([dfold,'fMRI',filesep,cname,filesep,'dataset_volsurf.mat']);
dinfo = loadjson([dfold,'fMRI',filesep,cname,filesep,'dataset_info.json']);
dataset.type = dinfo.type;
dataset.subsets = dinfo.subsets;

if exist('outputmask','var')
    if ~isempty(outputmask)
        omask=ea_load_nii(outputmask);
        omaskidx=find(omask.img(:));
    else
        omaskidx=dataset.vol.outidx;
    end
else
    omaskidx=dataset.vol.outidx; % use all.
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

        [~,seedfn{s,lr}]=fileparts(sfile{s,lr});
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
for s=1:numseed-1
    fX{s}=nan(length(omaskidx),numsub);
end

ea_dispercent(0,'Iterating through subjects');

scnt=1;
for mcfi=usesubjects % iterate across subjects
    howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi);
    targetix=sweightidx{1};
    clear stc
    thiscorr=cell(numseed-1,1);
    for s=1:numseed-1
        thiscorr{s}=zeros(length(omaskidx),howmanyruns);
    end
    for run=1:howmanyruns
        for s=2:numseed
            switch dataset.type
                case 'fMRI_matrix'
                    ea_error('Command partial map is not supported for matrix type datasets.')
                case 'fMRI_timecourses'
                    load([dfoldvol,dataset.vol.subIDs{mcfi}{run+1}])
                    gmtc=single(gmtc);
                    stc(:,s-1)=mean(gmtc(sweightidx{s},:).*repmat(sweightidxmx{s},1,size(gmtc,2)));
            end
        end

        % now we have all seeds, need to iterate across voxels of
        % target to get pmap values
        for s=1:size(stc,2)
            seedstc=stc(:,s);
            otherstc=stc;
            otherstc(:,s)=[];

            targtc=gmtc(targetix,:);
            thiscorr{s}(targetix,run)=partialcorr(targtc',seedstc,otherstc);

        end
    end

    for s=1:size(stc,2)
        fX{s}(:,scnt)=mean(thiscorr{s},2);
        if writeoutsinglefiles
            ccmap=dataset.vol.space;
            ccmap.dt=[16 0];
            ccmap.img=single(ccmap.img);

            if ~isempty(outputfolder)
                ccmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_corr.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    ccmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_id-',dataset.vol.subIDs{mcfi}{1},'_corr.nii']);
                else
                    ccmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'id', dataset.vol.subIDs{mcfi}{1}, 'suffix', 'corr');
                end
            end

            ccmap.img(omaskidx)=fX{s}(:,scnt);
            spm_write_vol(ccmap,ccmap.img);
        end
    end

    ea_dispercent(scnt/numsub);
    scnt=scnt+1;
end
ea_dispercent(1,'end');
seedfn(1)=[]; % delete first seed filename (which is target).

switch dataset.type
    case 'fMRI_matrix'
        ea_error('Command partial map is not supported for matrix type datasets.')
    case 'fMRI_timecourses'
        for s=1:size(seedfn,1) % subtract 1 in case of pmap command
            % export mean
            M=ea_nanmean(fX{s}',1);
            mmap=dataset.vol.space;
            mmap.dt(1) = 16;
            mmap.img(:)=0;
            mmap.img=single(mmap.img);
            mmap.img(omaskidx)=M;

            if ~isempty(outputfolder)
                mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-AvgR_funcpmap1stseed.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-AvgR_funcpmap1stseed.nii']);
                else
                    mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'AvgR', 'suffix', 'funcpmap1stseed');
                end
            end

            ea_write_nii(mmap);
            if usegzip
                gzip(mmap.fname);
                delete(mmap.fname);
            end

            % export variance
            M=ea_nanvar(fX{s}');
            mmap=dataset.vol.space;
            mmap.dt(1) = 16;
            mmap.img(:)=0;
            mmap.img=single(mmap.img);
            mmap.img(omaskidx)=M;

            if ~isempty(outputfolder)
                mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-VarR_funcpmap1stseed.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-VarR_funcpmap1stseed.nii']);
                else
                    mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'VarR', 'suffix', 'funcpmap1stseed');
                end
            end

            ea_write_nii(mmap);
            if usegzip
                gzip(mmap.fname);
                delete(mmap.fname);
            end

            % fisher-transform:
            fX{s}=atanh(fX{s});

            % export fz-mean
            M=nanmean(fX{s}');
            mmap=dataset.vol.space;
            mmap.dt(1) = 16;
            mmap.img(:)=0;
            mmap.img=single(mmap.img);
            mmap.img(omaskidx)=M;

            if ~isempty(outputfolder)
                mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-AvgRFz_funcpmap1stseed.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-AvgRFz_funcpmap1stseed.nii']);
                else
                    mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'AvgRFz', 'suffix', 'funcpmap1stseed');
                end
            end

            spm_write_vol(mmap,mmap.img);
            if usegzip
                gzip(mmap.fname);
                delete(mmap.fname);
            end

            % export T
            [~,~,~,tstat]=ttest(fX{s}');
            tmap=dataset.vol.space;
            tmap.img(:)=0;
            tmap.dt(1) = 16;
            tmap.img=single(tmap.img);

            tmap.img(omaskidx)=tstat.tstat;

            if ~isempty(outputfolder)
                tmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-T_funcpmap1stseed.nii']);
            else
                if ~isBIDSFileName(sfile{s})
                    tmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-T_funcpmap1stseed.nii']);
                else
                    tmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'T', 'suffix', 'funcpmap1stseed');
                end
            end

            spm_write_vol(tmap,tmap.img);
            if usegzip
                gzip(tmap.fname);
                delete(tmap.fname);
            end
        end
end

toc


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
