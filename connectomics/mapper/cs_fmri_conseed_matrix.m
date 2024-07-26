function cs_fmri_conseed_matrix(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,exportgmtc)

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

dfoldsurf=[dfold,'fMRI',filesep,cname,filesep,'surf',filesep];
dfoldvol=[dfold,'fMRI',filesep,cname,filesep,'vol',filesep]; % expand to /vol subdir.

dataset = load([dfold,'fMRI',filesep,cname,filesep,'dataset_volsurf.mat']);
dataset.connLabel = connLabel;
dinfo = loadjson([dfold,'fMRI',filesep,cname,filesep,'dataset_info.json']);
dataset.type = dinfo.type;
dataset.subsets = dinfo.subsets;

try
    dataset.formatversion=dinfo.formatversion;
catch
    dataset.formatversion=1.0;
end

if isempty(outputfolder)
    warning('off', 'backtrace');
    warning('Custom output folder not specified! Will save result to current folder.');
    warning('on', 'backtrace');
    outputfolder = [pwd, filesep];
elseif ~strcmp(outputfolder(end),filesep)
    outputfolder = [outputfolder,filesep];
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



% Check for Parallel Computing Toolbox
hasParallelToolbox = license('test', 'Distrib_Computing_Toolbox');

% Initialize vars
addp = '';
ea_dispercent(0, 'Iterating through subjects');
scnt = 1;

% Preallocate result arrays
fX = zeros(sum(sum(triu(ones(numseed), 1))), numsub);

% Temporary cell array for parfor loop
temp_fX = cell(1, numsub);

if hasParallelToolbox
    parfor mcfi_idx = 1:numel(usesubjects)
        mcfi = usesubjects(mcfi_idx);
        howmanyruns = ea_cs_dethowmanyruns(dataset, mcfi);

        thiscorr = zeros(sum(sum(triu(ones(numseed), 1))), howmanyruns);

        for run = 1:howmanyruns
            data = load_data(mcfi, dataset, dfoldvol, dfoldsurf, run, sfile);
            stc = compute_stc(data, sfile, sweightidx, sweightidxmx, numseed, dataset);

            if exportgmtc
                save_gmtc(stc, outputfolder, connLabel, dataset, mcfi, run);
            end

            thiscorr(:,run) = compute_correlation(stc, numseed, cmd, dataset);
        end

        thiscorr = mean(thiscorr, 2);
        temp_fX{mcfi_idx} = thiscorr;

        if writeoutsinglefiles
            X = zeros(numseed);
            X(:) = thiscorr;
            save_funcmatrix(X, outputfolder, connLabel, dataset, mcfi);
        end

        ea_dispercent(mcfi_idx / numel(usesubjects));
    end
else
    for mcfi_idx = 1:numel(usesubjects)
        mcfi = usesubjects(mcfi_idx);
        howmanyruns = ea_cs_dethowmanyruns(dataset, mcfi);

        thiscorr = zeros(sum(sum(triu(ones(numseed), 1))), howmanyruns);

        for run = 1:howmanyruns
            data = load_data(mcfi, dataset, dfoldvol, dfoldsurf, run, sfile);
            stc = compute_stc(data, sfile, sweightidx, sweightidxmx, numseed, dataset);

            if exportgmtc
                save_gmtc(stc, outputfolder, connLabel, dataset, mcfi, run);
            end

            thiscorr(:,run) = compute_correlation(stc, numseed, cmd, dataset);
        end

        thiscorr = mean(thiscorr, 2);
        fX(:, scnt) = thiscorr;

        if writeoutsinglefiles
            X = zeros(numseed);
            X(:) = thiscorr;
            save_funcmatrix(X, outputfolder, connLabel, dataset, mcfi);
        end

        ea_dispercent(scnt / numsub);
        scnt = scnt + 1;
    end
end

% Aggregate results from temporary cell array
if hasParallelToolbox
    for mcfi_idx = 1:numel(usesubjects)
        fX(:, mcfi_idx) = temp_fX{mcfi_idx};
    end
end

ea_dispercent(1, 'end');


switch dataset.type
    case 'fMRI_matrix'
        ea_error(['Command ',cmd,' in combination with an fMRI-matrix not (yet) supported.']);
    case 'fMRI_timecourses'
        seeds = sfile;
        if isfield(options.lcm, 'parcSeedName') && ~isempty(options.lcm.parcSeedName)
            isParcSeed = 1;
            parcSeedName = options.lcm.parcSeedName;
        else
            isParcSeed = 0;
        end

        % export mean
        M=nanmean(fX');
        X=zeros(numseed);
        X(logical(triu(ones(numseed),1)))=M;
        X=X+X';
        X(logical(eye(length(X))))=1;
        if isParcSeed
            save(fullfile(outputfolder, [parcSeedName, '_conn-', connLabel, '_desc-AvgR_func', cmd, '.mat']), 'X', '-v7.3');
        else
            if ~isBIDSFileName(sfile{1})
                [~, fn] = fileparts(sfile{1});
                save(fullfile(outputfolder, [fn, '_conn-', connLabel, '_desc-AvgR_func', cmd, '.mat']), 'X', 'seeds', '-v7.3');
            else
                save(setBIDSEntity(sfile{s}, 'dir', outputfolder, 'sub', '', 'conn', connLabel, 'desc', 'AvgR', 'suffix', ['func',cmd], 'ext', 'mat'), 'X', 'seeds','-v7.3');
            end
        end

        % export variance
        M=nanvar(fX');
        X=zeros(numseed);
        X(logical(triu(ones(numseed),1)))=M;
        X=X+X';
        X(logical(eye(length(X))))=1;
        if isParcSeed
            save(fullfile(outputfolder, [parcSeedName, '_conn-', connLabel, '_desc-VarR_func', cmd, '.mat']), 'X', '-v7.3');
        else
            if ~isBIDSFileName(sfile{1})
                [~, fn] = fileparts(sfile{1});
                save(fullfile(outputfolder, [fn, '_conn-', connLabel, '_desc-VarR_func', cmd, '.mat']), 'X', 'seeds', '-v7.3');
            else
                save(setBIDSEntity(sfile{s}, 'dir', outputfolder, 'sub', '', 'conn', connLabel, 'desc', 'VarR', 'suffix', ['func',cmd], 'ext', 'mat'), 'X', 'seeds','-v7.3');
            end
        end

        % fisher-transform:
        fX=atanh(fX);
        M=nanmean(fX');
        X=zeros(numseed);
        X(logical(triu(ones(numseed),1)))=M;
        X=X+X';
        X(logical(eye(length(X))))=1;
        if isParcSeed
            save(fullfile(outputfolder, [parcSeedName, '_conn-', connLabel, '_desc-AvgRFz_func', cmd, '.mat']), 'X', '-v7.3');
        else
            if ~isBIDSFileName(sfile{1})
                [~, fn] = fileparts(sfile{1});
                save(fullfile(outputfolder, [fn, '_conn-', connLabel, '_desc-AvgRFz_func', cmd, '.mat']), 'X', 'seeds', '-v7.3');
            else
                save(setBIDSEntity(sfile{s}, 'dir', outputfolder, 'sub', '', 'conn', connLabel, 'desc', 'AvgRFz', 'suffix', ['func',cmd], 'ext', 'mat'), 'X', 'seeds','-v7.3');
            end
        end

        % export T
        [~,~,~,tstat]=ttest(fX');
        X=zeros(numseed);
        X(logical(triu(ones(numseed),1)))=tstat.tstat;
        X=X+X';
        X(logical(eye(length(X))))=1;
        if isParcSeed
            save(fullfile(outputfolder, [parcSeedName, '_conn-', connLabel, '_desc-T_func', cmd, '.mat']), 'X', '-v7.3');
        else
            if ~isBIDSFileName(sfile{1})
                [~, fn] = fileparts(sfile{1});
                save(fullfile(outputfolder, [fn, '_conn-', connLabel, '_desc-T_func', cmd, '.mat']), 'X', 'seeds', '-v7.3');
            else
                save(setBIDSEntity(sfile{s}, 'dir', outputfolder, 'sub', '', 'conn', connLabel, 'desc', 'T', 'suffix', ['func',cmd], 'ext', 'mat'), 'X', 'seeds','-v7.3');
            end
        end
end
toc



function data = load_data(mcfi, dataset, dfoldvol, dfoldsurf, run, sfile)
data = struct();
data.gmtc = load([dfoldvol, dataset.vol.subIDs{mcfi}{run+1}], 'gmtc').gmtc;
data.gmtc = single(data.gmtc);

if size(sfile, 2) > 1
    data.ls = load([dfoldsurf, dataset.surf.l.subIDs{mcfi}{run+1}]);
    data.rs = load([dfoldsurf, dataset.surf.r.subIDs{mcfi}{run+1}]);
    data.ls.gmtc = single(data.ls.gmtc);
    data.rs.gmtc = single(data.rs.gmtc);
end





function stc = compute_stc(data, sfile, sweightidx, sweightidxmx, numseed, dataset)
if dataset.formatversion == 1
    stc = nan(numseed, size(data.gmtc, 2));
elseif dataset.formatversion > 1
    stc = nan(size(data.gmtc, 1), numseed);
end

for s = 1:numseed
    if size(sfile, 2) > 1
        stc(s,:) = mean([data.ls.gmtc(sweightidx{s,1},:) .* repmat(sweightidxmx{s,1}, 1, size(data.ls.gmtc, 2)); ...
            data.rs.gmtc(sweightidx{s,2},:) .* repmat(sweightidxmx{s,2}, 1, size(data.rs.gmtc, 2))], 1);
    else
        if dataset.formatversion == 1
            stc(s,:) = mean(data.gmtc(sweightidx{s}, :) .* repmat(sweightidxmx{s}, 1, size(data.gmtc, 2)), 1);
        elseif dataset.formatversion > 1
            stc(:,s) = mean(data.gmtc(:, sweightidx{s}) .* repmat(sweightidxmx{s}', size(data.gmtc, 1), 1), 2)';
        end
    end
end




function thiscorr = compute_correlation(stc, numseed, cmd, dataset)
if strcmp(cmd, 'matrix')
    if dataset.formatversion == 1
        X = corrcoef(stc');
    elseif dataset.formatversion > 1
        X = corrcoef(stc);
    end
elseif strcmp(cmd, 'pmatrix')
    if dataset.formatversion == 1
        X = partialcorr(stc');
    elseif dataset.formatversion > 1
        X = partialcorr(stc);
    end
end
thiscorr = X(logical(triu(ones(numseed), 1)));




function save_gmtc(stc, outputfolder, connLabel, dataset, mcfi, run)
tmp.gmtc = stc;
save(fullfile(outputfolder, ['conn-', connLabel, '_id-', dataset.vol.subIDs{mcfi}{1}, '_run-', num2str(run, '%02d'), '_timeseries.mat']), '-struct', 'tmp', '-v7.3');


function save_funcmatrix(X, outputfolder, connLabel, dataset, mcfi)
save(fullfile(outputfolder, ['conn-', connLabel, '_id-', dataset.vol.subIDs{mcfi}{1}, '_funcmatrix.mat']), 'X', '-v7.3');


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
