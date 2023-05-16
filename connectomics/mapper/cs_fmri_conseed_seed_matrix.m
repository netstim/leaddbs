function cs_fmri_conseed_seed_matrix(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask)

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

dataset = load([dfold,'fMRI',filesep,cname,filesep,'dataset_volsurf.mat']);
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
        sweightidxmx{s,lr}=sweightidxmx{s,lr}./sum(sweightidxmx{s,lr});
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

numSubUse=length(dataset.vol.subIDs);

if ~exist('subset','var') % use all subjects
    usesubjects=1:numSubUse;
else
    for ds=1:length(dataset.subsets)
        if strcmp(subset,dataset.subsets{ds}.name)
            usesubjects=dataset.subsets{ds}.subs;
            break
        end
    end
    numSubUse=length(usesubjects);
end

if ~exist('db','var')
    db=matfile([dfold,'fMRI',filesep,cname,filesep,'AllX.mat'],'Writable',false);
    probe=db.X(1,1);
    switch class(probe)
        case 'int16'
            needdivide=1;
        case {'single','double'}
            needdivide=0;
        otherwise
            ea_error('File format not supported');
    end
end

stack = dbstack;
if ismember('ea_networkmapping_calcvals', {stack.name})
    isNetworkMappingRun = 1;
else
    isNetworkMappingRun = 0;
end

disp('Iterating through subjects...');
for subj = 1:numSubUse % iterate across subjects
    mcfi = usesubjects(subj);
    howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi);

    for s=1:numseed
        for run=1:howmanyruns
            Rw=nan(length(sweightidx{s}),pixdim);
            ea_dispercent(0,'Parsing connectome');
            queryfrom=1;
            while 1 % only possible via loop given matfile mapping restrictions, querying as efficiently as possible.
                queryuntil=queryfrom;
                for queryinterval=1:length(sweightidx{s})-queryfrom
                    if ~isequal(sweightidx{s}(queryfrom)+queryinterval,sweightidx{s}(queryfrom+queryinterval)) % check if continuous
                       break
                    end
                    queryuntil=queryfrom+queryinterval;
                end
                Rw(queryfrom:queryuntil,:)=db.X(sweightidx{s}(queryfrom:queryuntil),1:pixdim);
                ea_dispercent(queryuntil/length(sweightidx{s}));
                queryfrom=queryuntil+1;
                if queryfrom>length(sweightidx{s})
                    break
                end
            end
            ea_dispercent(1,'end');
            if needdivide
                Rw=Rw/(2^15); % convert to actual R values
            end
            Rw=Rw.*repmat(sweightidxmx{s},1,pixdim); % map weights of seed to entries
            Rw=sum(Rw,1); % sum is fine since sum of sweightidxmx{s} == 1
        end

        mmap=dataset.vol.space;

        if ~isempty(outputfolder)
            mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-AvgR_funcmap.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-AvgR_funcmap.nii']);
            else
                mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'AvgR', 'suffix', 'funcmap');
            end
        end

        mmap.dt(1) = 16;
        mmap.img(:)=0;
        mmap.img=single(mmap.img);
        mmap.img(omaskidx)=Rw;

        if ~isNetworkMappingRun
            ea_write_nii(mmap);
            if usegzip
                gzip(mmap.fname);
                delete(mmap.fname);
            end
        end

        if ~isempty(outputfolder)
            mmap.fname = fullfile(outputfolder, [seedfn{s}, '_conn-', connLabel, '_desc-AvgRFz_funcmap.nii']);
        else
            if ~isBIDSFileName(sfile{s})
                mmap.fname = strrep(sfile{s}, '.nii', ['_conn-', connLabel, '_desc-AvgRFz_funcmap.nii']);
            else
                mmap.fname = setBIDSEntity(sfile{s}, 'conn', connLabel, 'desc', 'AvgRFz', 'suffix', 'funcmap');
            end
        end

        mmap.img(:)=0;
        mmap.img=single(mmap.img);
        mmap.img(omaskidx)=atanh(Rw);
        ea_write_nii(mmap);
        if usegzip
            gzip(mmap.fname);
            delete(mmap.fname);
        end
    end
end
disp('Done.');

toc


function howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi)
if strcmp(dataset.type,'fMRI_matrix')
    howmanyruns=1;
else
    howmanyruns=length(dataset.vol.subIDs{mcfi})-1;
end
