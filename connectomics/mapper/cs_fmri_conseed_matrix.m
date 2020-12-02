function cs_fmri_conseed_matrix(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,exportgmtc)

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
ocname=cname;

if ismember('>',cname)
    delim=strfind(cname,'>');
    subset=cname(delim+1:end);
    cname=cname(1:delim-1);
end

prefs = ea_prefs;
dfoldsurf=[dfold,'fMRI',filesep,cname,filesep,'surf',filesep];
dfoldvol=[dfold,'fMRI',filesep,cname,filesep,'vol',filesep]; % expand to /vol subdir.

d=load([dfold,'fMRI',filesep,cname,filesep,'dataset_info.mat']);
dataset=d.dataset;
clear d;
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

owasempty=0;
if ~exist('outputfolder','var')
    outputfolder=ea_getoutputfolder(sfile,ocname);
    owasempty=1;
else
    if isempty(outputfolder) % from shell wrapper.
        outputfolder=ea_getoutputfolder(sfile,ocname);
        owasempty=1;
    end
    if ~strcmp(outputfolder(end),filesep)
        outputfolder=[outputfolder,filesep];
    end
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
        [ixx]=unique(round(sweights)); ixx(ixx==0)=[];
        numseed=length(ixx);
        for parcseed=ixx'
            sweightidx{parcseed+1,1}=sweightidx{1,1}(sweightidxmx{1,1}==parcseed);
            sweightidxmx{parcseed+1,1}=ones(size(sweightidx{parcseed+1,1},1),1);
        end
        sweightidx(1)=[]; % original parcellation which has now been expanded to single seeds
        sweightidxmx(1)=[]; % original parcellation which has now been expanded to single seeds
        sfile=repmat(sfile,size(sweightidx,1),1);
    end
end

disp([num2str(numseed),' seeds, command = ',cmd,'.']);

numSubUse = length(dataset.vol.subIDs);

if ~exist('subset','var') % use all subjects
    usesubjects = 1:numSubUse;
else
    for ds=1:length(dataset.subsets)
        if strcmp(subset,dataset.subsets(ds).name)
            usesubjects=dataset.subsets(ds).subs;
            break
        end
    end
    numSubUse = length(usesubjects);
end

fX=nan(((numseed^2)-numseed)/2,numSubUse);

% init vars:
addp='';

isSurfAvail = isfield(dataset,'surf');
includeSurf = prefs.lcm.includesurf;

disp('Iterating through subjects...');
parfor subj = 1:numSubUse % iterate across subjects
    mcfi = usesubjects(subj);
    disp(['Subject ', num2str(mcfi, '%04d'),'/',num2str(numSubUse,'%04d'),'...']);
    howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi);

    subIDVol = dataset.vol.subIDs{mcfi};
    if isSurfAvail && includeSurf
        subIDSurfL = dataset.surf.l.subIDs{mcfi};
        subIDSurfR = dataset.surf.r.subIDs{mcfi};
    end

    thiscorr = nan(numseed^2,howmanyruns);

    for run=1:howmanyruns
        gmtcstruc = load([dfoldvol,subIDVol{run+1}]);
        gmtc = double(gmtcstruc.gmtc);

        if isSurfAvail && includeSurf
            ls_struc = load([dfoldsurf,subIDSurfL{run+1}]);
            rs_struc = load([dfoldsurf,subIDSurfR{run+1}]);
            ls_gmtc = double(ls_struc.gmtc);
            rs_gmtc = double(rs_struc.gmtc);
        end

        stc = nan(numseed,size(gmtc,2));

        for s=1:numseed
            if size(sfile(s,:),2)>1 % dealing with surface seed
                stc(s,:) = mean([ls_gmtc(sweightidx{s,1},:) .* repmat(sweightidxmx{s,1},1,size(ls_gmtc,2));...
                                 rs_gmtc(sweightidx{s,2},:) .* repmat(sweightidxmx{s,2},1,size(rs_gmtc,2))] ,1); % seed time course
            else % volume seed
                try
                    stc(s,:) = mean(gmtc(sweightidx{s},:) .* repmat(sweightidxmx{s},1,size(gmtc,2)), 1); % seed time course
                catch
                    keyboard
                end
            end
        end

        if exportgmtc
            mat = matfile([outputfolder,addp,'gmtc_',subIDVol{1},'_run',num2str(run,'%02d'),'.mat'], 'Writable', true);
            mat.gmtc = stc;
        end

        switch cmd
            case 'matrix'
                X = corrcoef(stc');
            case 'pmatrix'
                X = partialcorr(stc');
        end
        thiscorr(:,run) = X(:);

    end

    thiscorrAvg = mean(thiscorr,2);
    X(:) = thiscorrAvg;
    fX(:, subj) = X(logical(triu(ones(numseed),1)));

    if writeoutsinglefiles
        mat = matfile([outputfolder,addp,'corrMx_',subIDVol{1},'.mat'], 'Writable', true);
        mat.X = X;
    end
end
disp('Done.');

% export mean
M=nanmean(fX');
X=zeros(numseed);
X(logical(triu(ones(numseed),1)))=M;
X=X+X';
X(logical(eye(length(X))))=1;
save([outputfolder,cmd,'_corrMx_AvgR.mat'],'X','-v7.3');

% export variance
M=nanvar(fX');
X=zeros(numseed);
X(logical(triu(ones(numseed),1)))=M;
X=X+X';
X(logical(eye(length(X))))=1;
save([outputfolder,cmd,'_corrMx_VarR.mat'],'X','-v7.3');

% fisher-transform:
fX=atanh(fX);
M=nanmean(fX');
X=zeros(numseed);
X(logical(triu(ones(numseed),1)))=M;
X=X+X';
X(logical(eye(length(X))))=1;
save([outputfolder,cmd,'_corrMx_AvgR_Fz.mat'],'X','-v7.3');

% export T
[~,~,~,tstat]=ttest(fX');
X=zeros(numseed);
X(logical(triu(ones(numseed),1)))=tstat.tstat;
X=X+X';
X(logical(eye(length(X))))=1;
save([outputfolder,cmd,'_corrMx_T.mat'],'X','-v7.3');

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
if rightmat==loaded;
    return
end

load([datadir,num2str(rightmat),'.mat']);
loaded=rightmat;
