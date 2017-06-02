function cs_fmri_conseed(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask)
tic


if ~isdeployed
    addpath(genpath('/autofs/cluster/nimlab/connectomes/software/lead_dbs'));
    addpath('/autofs/cluster/nimlab/connectomes/software/spm12');
end
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

if ismember('>',cname)
    delim=strfind(cname,'>');
    subset=cname(delim+1:end);
    cname=cname(1:delim-1);
end


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


[sfile,roilist]=ea_handleseeds(sfile);



if ~exist('outputfolder','var')
    ea_getoutputfolder(sfile)
else
    if isempty(outputfolder) % from shell wrapper.
        ea_getoutputfolder(sfile)
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

for s=1:length(sfile)
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
        % assure sum of sweights is 1
        %sweights(logical(sweights))=sweights(logical(sweights))/abs(sum(sweights(logical(sweights))));
        sweightmx=repmat(sweights,1,1);
        
        sweightidx{s,lr}=find(sweights);
        sweightidxmx{s,lr}=double(sweightmx(sweightidx{s,lr},:));
    end
end

numseed=s;
disp([num2str(numseed),' seeds, command = ',cmd,'.']);


pixdim=length(dataset.vol.outidx);

numsub=length(dataset.vol.subIDs);
% init vars:
switch cmd
    case {'seed','pseed'}
        for s=1:numseed
            fX{s}=nan(length(omaskidx),numsub);
            rh.fX{s}=nan(10242,numsub);
            lh.fX{s}=nan(10242,numsub);
        end

    case 'pmap'
        
        for s=1:numseed-1
            fX{s}=nan(length(omaskidx),numsub);
        end
    otherwise
        fX=nan(((numseed^2)-numseed)/2,numsub);
end

switch cmd
    case 'matrix'
        addp='';
    case 'pmatrix'
        addp='p';
end



ea_dispercent(0,'Iterating through subjects');

if ~exist('subset','var') % use all subjects
    usesubjects=1:numsub;
else
    for ds=1:length(dataset.subsets)
        if strcmp(subset,dataset.subsets(ds).name)
            usesubjects=dataset.subsets(ds).subs;
            break
        end
    end
end



for mcfi=usesubjects % iterate across subjects
    howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi);
    switch cmd
        
        case {'seed','pseed'}
            
            for s=1:numseed
                
                thiscorr=zeros(length(omaskidx),howmanyruns);
                
                for run=1:howmanyruns
                    switch dataset.type
                        case 'fMRI_matrix'
                            if strcmp(cmd,'pseed')
                                ea_error('Cannot run partial seed on fMRI Matrix dataset.');
                            end

                            Rw=nan(length(sweightidx{s}),pixdim);
                            
                            if ~exist('db','var')
                                db=matfile([dfold,'fMRI',filesep,cname,filesep,'AllX.mat']);
                            end
                            
                           
                            cnt=1;
                            for ix=sweightidx{s}'
                                %    testnii.img(outidx)=mat(entry,:); % R
                                Rw(cnt,:)=db.X(sweightidx{s}(cnt),:);
                                cnt=cnt+1;
                            end
                            Rw=mean(Rw,1);
                            Rw=Rw/(2^15);
                        case 'fMRI_timecourses'
                            load([dfoldvol,dataset.vol.subIDs{mcfi}{run+1}])
                            if isfield(dataset,'surf')
                                % include surface:
                                ls=load([dfoldsurf,dataset.surf.l.subIDs{mcfi}{run+1}]);
                                rs=load([dfoldsurf,dataset.surf.r.subIDs{mcfi}{run+1}]);
                            end
                            gmtc=single(gmtc);
                            
                            switch cmd % build up seed tc for present subject
                                case 'seed'
                                    if size(sfile(s,:),2)>1 % dealing with surface seed
                                        ls.gmtc=single(ls.gmtc); rs.gmtc=single(rs.gmtc);
                                        stc=mean([ls.gmtc(sweightidx{s,1},:).*repmat(sweightidxmx{s,1},1,size(ls.gmtc,2));...
                                            rs.gmtc(sweightidx{s,2},:).*repmat(sweightidxmx{s,2},1,size(ls.gmtc,2))],1); % seed time course
                                    else % volume seed

                                        stc=mean(gmtc(sweightidx{s},:).*repmat(sweightidxmx{s},1,size(gmtc,2)),1); % seed time course

                                    end
                                    
                                    
                                case 'pseed'
                                    clear stc
                                    for subseed=1:numseed
                                        if size(sfile(subseed,:),2)>1 % dealing with surface seed
                                            ls.gmtc=single(ls.gmtc); rs.gmtc=single(rs.gmtc);
                                            stc(:,subseed)=mean([ls.gmtc(sweightidx{subseed,1},:).*repmat(sweightidxmx{subseed,1},1,size(ls.gmtc,2));...
                                                rs.gmtc(sweightidx{subseed,2},:).*repmat(sweightidxmx{subseed,2},1,size(rs.gmtc,2))],1); % seed time course
                                        else % volume seed
                                            stc(:,subseed)=mean(gmtc(sweightidx{subseed},:).*repmat(sweightidxmx{subseed},1,size(gmtc,2)),1); % seed time course
                                        end
                                    end
                                    os=1:numseed; os(s)=[]; % remaining seeds
                                    [~,~,stc]=regress(stc(:,s),addone(stc(:,os))); % regress out other time series from current one
                                    stc=stc';
                            end
                            thiscorr(:,run)=corr(stc',gmtc(maskuseidx,:)','type','Pearson');
                            if isfield(dataset,'surf')
                                % include surface:
                                rs.gmtc=single(rs.gmtc);
                                ls.gmtc=single(ls.gmtc);
                                ls.thiscorr(:,run)=corr(stc',ls.gmtc','type','Pearson');
                                rs.thiscorr(:,run)=corr(stc',rs.gmtc','type','Pearson');
                            end
                    end
                end
                
                fX{s}(:,mcfi)=mean(thiscorr,2);
                if isfield(dataset,'surf')
                    lh.fX{s}(:,mcfi)=mean(ls.thiscorr,2);
                    rh.fX{s}(:,mcfi)=mean(rs.thiscorr,2);
                end
                
                if writeoutsinglefiles && (~strcmp(dataset.type,'fMRI_matrix'))
                    ccmap=dataset.vol.space;
                    ccmap.img=single(ccmap.img);
                    ccmap.fname=[outputfolder,seedfn{s},'_',dataset.vol.subIDs{mcfi}{1},'_corr.nii'];
                    ccmap.img(omaskidx)=fX{s}(:,mcfi);
                    ccmap.dt=[16,0];
                    spm_write_vol(ccmap,ccmap.img);
                end
            end

            
            
        case 'pmap'
            
            
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
                            
                            ea_error('Pmap not supported with use of fMRI_Matrix (yet).');
                            
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
                fX{s}(:,mcfi)=mean(thiscorr{s},2);
                if writeoutsinglefiles
                    ccmap=dataset.vol.space;
                    ccmap.dt=[16 0];
                    ccmap.img=single(ccmap.img);
                    ccmap.fname=[outputfolder,seedfn{s},'_',dataset.vol.subIDs{mcfi}{1},'_pmap.nii'];
                    ccmap.img(omaskidx)=fX{s}(:,mcfi);
                    spm_write_vol(ccmap,ccmap.img);
                end
            end
            
            
        otherwise
            clear stc
            for run=1:howmanyruns
                load([dfoldvol,dataset.vol.subIDs{mcfi}{run+1}])
                gmtc=single(gmtc);
                
                if size(sfile(s,:),2)>1
                    % include surface:
                    ls=load([dfoldsurf,dataset.surf.l.subIDs{mcfi}{run+1}]);
                    rs=load([dfoldsurf,dataset.surf.r.subIDs{mcfi}{run+1}]);
                    ls.gmtc=single(ls.gmtc); rs.gmtc=single(rs.gmtc);
                end
                for s=1:numseed
                    if size(sfile(s,:),2)>1 % dealing with surface seed
                        stc(s,:)=mean([ls.gmtc(sweightidx{s,1},:).*repmat(sweightidxmx{s,1},1,size(ls.gmtc,2));...
                            rs.gmtc(sweightidx{s,2},:).*repmat(sweightidxmx{s,2},1,size(rs.gmtc,2))],1); % seed time course
                    else % volume seed
                        try
                        stc(s,:)=mean(gmtc(sweightidx{s},:).*repmat(sweightidxmx{s},1,size(gmtc,2)),1); % seed time course
                        catch
                            keyboard
                        end
                    end
                end
                
                switch cmd
                    case 'matrix'
                        X=corrcoef(stc');
                        
                    case 'pmatrix'
                        X=partialcorr(stc');
                end
                thiscorr(:,run)=X(:);
                
            end
            thiscorr=mean(thiscorr,2);
            X(:)=thiscorr;
            fX(:,mcfi)=X(logical(triu(ones(numseed),1)));
            
            if writeoutsinglefiles
                save([outputfolder,addp,'corrMx_',dataset.vol.subIDs{mcfi}{1},'.mat'],'X','-v7.3');
            end
    end
        ea_dispercent(mcfi/numsub);
end
ea_dispercent(1,'end');
ispmap=strcmp(cmd,'pmap');
if ispmap
    seedfn(1)=[]; % delete first seed filename (which is target).
end
switch dataset.type
    case 'fMRI_matrix'
        switch cmd
            case {'seed'}
                
                mmap=dataset.vol.space;
                mmap.dt=[16,0];
                mmap.img(:)=0;
                mmap.img=single(mmap.img);
                mmap.img(omaskidx)=Rw;
                mmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR.nii'];
                ea_write_nii(mmap);
                if usegzip
                    gzip(mmap.fname);
                    delete(mmap.fname);
                end
                
                
                
            otherwise
                ea_error(['Command ',cmd,' in combination with an fMRI-matrix not (yet) supported.']);
        end
        
    case 'fMRI_timecourses'
        switch cmd
            case {'seed','pmap','pseed'}
                for s=1:length(seedfn) % subtract 1 in case of pmap command
                    
                    % export mean
                    M=nanmean(fX{s}');
                    mmap=dataset.vol.space;
                    mmap.dt=[16,0];
                    mmap.img(:)=0;
                    mmap.img=single(mmap.img);
                    mmap.img(omaskidx)=M;
                    
                    mmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR.nii'];
                    ea_write_nii(mmap);
                    if usegzip
                        gzip(mmap.fname);
                        delete(mmap.fname);
                    end
                    
                    % export variance
                    M=nanvar(fX{s}');
                    mmap=dataset.vol.space;
                    mmap.dt=[16,0];
                    mmap.img(:)=0;
                    mmap.img=single(mmap.img);
                    mmap.img(omaskidx)=M;
                    
                    mmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_VarR.nii'];
                    ea_write_nii(mmap);
                    if usegzip
                        gzip(mmap.fname);
                        delete(mmap.fname);
                    end
                    
                    if ~ispmap && isfield(dataset,'surf')
                        % lh surf
                        lM=nanmean(lh.fX{s}');
                        lmmap=dataset.surf.l.space;
                        lmmap.dt=[16,0];
                        lmmap.img=zeros([size(lmmap.img,1),size(lmmap.img,2),size(lmmap.img,3)]);
                        lmmap.img=single(lmmap.img);
                        lmmap.img(:)=lM(:);
                        lmmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR_surf_lh.nii'];
                        ea_write_nii(lmmap);
                        if usegzip
                            gzip(lmmap.fname);
                            delete(lmmap.fname);
                        end
                        
                        % rh surf
                        rM=nanmean(rh.fX{s}');
                        rmmap=dataset.surf.r.space;
                        rmmap.dt=[16,0];
                        rmmap.img=zeros([size(rmmap.img,1),size(rmmap.img,2),size(rmmap.img,3)]);
                        rmmap.img=single(rmmap.img);
                        rmmap.img(:)=rM(:);
                        rmmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR_surf_rh.nii'];
                        ea_write_nii(rmmap);
                        if usegzip
                            gzip(rmmap.fname);
                            delete(rmmap.fname);
                        end
                    end
                    
                    
                    % fisher-transform:
                    fX{s}=atanh(fX{s});
                    if ~ispmap && isfield(dataset,'surf')
                        lh.fX{s}=atanh(lh.fX{s});
                        rh.fX{s}=atanh(rh.fX{s});
                    end
                    % export fz-mean
                    
                    M=nanmean(fX{s}');
                    mmap=dataset.vol.space;
                    mmap.dt=[16,0];
                    mmap.img(:)=0;
                    mmap.img=single(mmap.img);
                    mmap.img(omaskidx)=M;
                    mmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR_Fz.nii'];
                    spm_write_vol(mmap,mmap.img);
                    if usegzip
                        gzip(mmap.fname);
                        delete(mmap.fname);
                    end
                    if ~ispmap && isfield(dataset,'surf')
                        % lh surf
                        lM=nanmean(lh.fX{s}');
                        lmmap=dataset.surf.l.space;
                        lmmap.dt=[16,0];
                        lmmap.img=zeros([size(lmmap.img,1),size(lmmap.img,2),size(lmmap.img,3)]);
                        lmmap.img=single(lmmap.img);
                        lmmap.img(:)=lM(:);
                        lmmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR_Fz_surf_lh.nii'];
                        ea_write_nii(lmmap);
                        if usegzip
                            gzip(lmmap.fname);
                            delete(lmmap.fname);
                        end
                        
                        % rh surf
                        rM=nanmean(rh.fX{s}');
                        rmmap=dataset.surf.r.space;
                        rmmap.dt=[16,0];
                        rmmap.img=zeros([size(rmmap.img,1),size(rmmap.img,2),size(rmmap.img,3)]);
                        rmmap.img=single(rmmap.img);
                        rmmap.img(:)=rM(:);
                        rmmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_AvgR_Fz_surf_rh.nii'];
                        ea_write_nii(rmmap);
                        if usegzip
                            gzip(rmmap.fname);
                            delete(rmmap.fname);
                        end
                    end
                    
                    % export T
                    
                    [~,~,~,tstat]=ttest(fX{s}');
                    tmap=dataset.vol.space;
                    tmap.img(:)=0;
                    tmap.dt=[16,0];
                    tmap.img=single(tmap.img);
                    
                    tmap.img(omaskidx)=tstat.tstat;
                    
                    tmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_T.nii'];
                    spm_write_vol(tmap,tmap.img);
                    if usegzip
                        gzip(tmap.fname);
                        delete(tmap.fname);
                    end
                    
                    
                    
                    
                    
                    if ~ispmap && isfield(dataset,'surf')
                        % lh surf
                        [~,~,~,ltstat]=ttest(lh.fX{s}');
                        lmmap=dataset.surf.l.space;
                        lmmap.dt=[16,0];
                        lmmap.img=zeros([size(lmmap.img,1),size(lmmap.img,2),size(lmmap.img,3)]);
                        lmmap.img=single(lmmap.img);
                        lmmap.img(:)=ltstat.tstat(:);
                        lmmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_T_surf_lh.nii'];
                        ea_write_nii(lmmap);
                        if usegzip
                            gzip(lmmap.fname);
                            delete(lmmap.fname);
                        end
                        
                        % rh surf
                        [~,~,~,rtstat]=ttest(rh.fX{s}');
                        rmmap=dataset.surf.r.space;
                        rmmap.dt=[16,0];
                        rmmap.img=zeros([size(rmmap.img,1),size(rmmap.img,2),size(rmmap.img,3)]);
                        rmmap.img=single(rmmap.img);
                        rmmap.img(:)=rtstat.tstat(:);
                        rmmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_T_surf_rh.nii'];
                        ea_write_nii(rmmap);
                        if usegzip
                            gzip(rmmap.fname);
                            delete(rmmap.fname);
                        end
                    end
                end
                
            otherwise
                
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
                
        end
end


toc


function s=ea_conformseedtofmri(dataset,s)
td=tempdir;
dataset.vol.space.fname=[td,'tmpspace.nii'];
ea_write_nii(dataset.vol.space);
s.fname=[td,'tmpseed.nii'];
ea_write_nii(s);

ea_conformspaceto([td,'tmpspace.nii'],[td,'tmpseed.nii']);
s=ea_load_nii(s.fname);


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

function outputfolder=ea_getoutputfolder(sfile)
file=sfile{1};
[pth,fn,ext]=fileparts(file); % exit to same folder as seed.
outputfolder=[pth,filesep];