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

dfoldsurf=[dfold,'fMRI',filesep,cname,filesep,'surf',filesep];
dfoldvol=[dfold,'fMRI',filesep,cname,filesep,'vol',filesep]; % expand to /vol subdir.



msk=ea_load_nii([dfold,'spacedefinitions',filesep,'222.nii']);
lmsk=ea_load_nii([dfold,'spacedefinitions',filesep,'lh_fsaverage.nii']);
rmsk=ea_load_nii([dfold,'spacedefinitions',filesep,'rh_fsaverage.nii']);

load([dfoldvol,'outidx.mat']);
load([dfoldvol,'subIDs.mat']);
load([dfoldsurf,'subIDs.mat']);
if exist('outputmask','var')
    if ~isempty(outputmask)
    omask=ea_load_nii(outputmask);
    omaskidx=find(omask.img(:));
    [~,maskuseidx]=ismember(omaskidx,outidx);
    else
        omaskidx=outidx;
        maskuseidx=1:length(outidx);
    end
else
    omaskidx=outidx; % use all.
        maskuseidx=1:length(outidx);
end

if iscell(sfile) % already supplied in cell format
    if length(sfile)>1
        roilist=1;
    else
        roilist=0;
    end
    
else
    [pth,fn,ext]=fileparts(sfile);
    if strcmp(ext,'.txt')
        roilist=1;
        
        sfile=getrois(sfile);
    else
        roilist=0;
        sfile={sfile};
    end
end



if ~exist('outputfolder','var')
    [pth,fn,ext]=fileparts(sfile); % exit to same folder as seed.
    outputfolder=[pth,filesep];
else
    if isempty(outputfolder) % from shell wrapper.
    [pth,fn,ext]=fileparts(sfile); % exit to same folder as seed.
    outputfolder=[pth,filesep];    
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
    seed{s}=ea_load_nii(ea_niigz(sfile{s}));
    
    [~,seedfn{s}]=fileparts(sfile{s});
    
    sweights=seed{s}.img(outidx);
    sweights(isnan(sweights))=0;
    sweights=double(sweights);
    % assure sum of sweights is 1
    sweights(logical(sweights))=sweights(logical(sweights))/sum(sweights(logical(sweights)));
    sweightmx=repmat(sweights,1,120);
    
    sweightidx{s}=find(sweights);
    sweightidxmx{s}=double(sweightmx(sweightidx{s},:));
end
numseed=s;


pixdim=length(outidx);

numsub=length(subIDs);
switch cmd
    case {'seed','seedvox_ram','seedvox_noram'}
        for s=1:numseed
            fX{s}=nan(length(omaskidx),numsub);
            rh.fX{s}=nan(10242,numsub);
            lh.fX{s}=nan(10242,numsub);
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


disp([num2str(numseed),' seeds, command = ',cmd,'.']);

ea_dispercent(0,'Iterating through subjects');

for mcfi=1:numsub
    ea_dispercent(mcfi/numsub);
    howmanyruns=length(subIDs{mcfi})-1;
    switch cmd
        
        case 'seed'
            
            for s=1:numseed
                thiscorr=zeros(length(omaskidx),howmanyruns);
                for run=1:howmanyruns
                    load([dfoldvol,subIDs{mcfi}{run+1}])
                    gmtc=single(gmtc);
                    stc=mean(gmtc(sweightidx{s},:).*sweightidxmx{s});
                    thiscorr(:,run)=corr(stc',gmtc(maskuseidx,:)','type','Pearson');
                    
                    % include surface:
                    ls=load([dfoldsurf,lsurfIDs{mcfi}{run+1}]);
                    rs=load([dfoldsurf,rsurfIDs{mcfi}{run+1}]);
                    rs.gmtc=single(rs.gmtc);
                    ls.gmtc=single(ls.gmtc);
                    ls.thiscorr(:,run)=corr(stc',ls.gmtc','type','Pearson');
                    rs.thiscorr(:,run)=corr(stc',rs.gmtc','type','Pearson');
                end
                
                fX{s}(:,mcfi)=mean(thiscorr,2);
                lh.fX{s}(:,mcfi)=mean(ls.thiscorr,2);
                rh.fX{s}(:,mcfi)=mean(rs.thiscorr,2);                
                
                
                if writeoutsinglefiles
                    ccmap=msk;
                    ccmap.img=single(ccmap.img);
                    ccmap.fname=[outputfolder,seedfn{s},'_',subIDs{mcfi}{1},'_corr.nii'];
                    ccmap.img(omaskidx)=fX{s}(:,mcfi);    
                    spm_write_vol(ccmap,ccmap.img);
                end
            end
            
        case {'seedvox_ram','seedvox_noram'}
            
               for s=1:numseed
                   swlength=length(sweightidx{s});

                       thiscorr=zeros(length(omaskidx),howmanyruns);
                for run=1:howmanyruns
                    load([dfoldvol,subIDs{mcfi}{run+1}])
                    gmtc=single(gmtc);
                    
                    switch cmd
                        case 'seedvox_ram'
                            stc=gmtc(logical(sweights),:);
                            X=corr(stc',gmtc');
                            Xweights=repmat(sweights(logical(sweights)),1,pixdim);
                            X=X.*Xweights;
                            clear Xweights
                            thiscorr(:,run)=mean(X,1);
                            clear X
                        case 'seedvox_noram'
                            gmtc=gmtc';
                            thiscorr=thiscorr';
                            
                            for seedvox=1:swlength
                                stc=gmtc(:,sweightidx{s}(seedvox));
                                addval=(corr(stc,gmtc(:,maskuseidx))*sweightidxmx{s}(seedvox,1));
                                addval(isnan(addval))=0;
                                thiscorr(run,:)=thiscorr(run,:)+addval;
                            end
                            
                            thiscorr=thiscorr';
                           thiscorr(:,run)=thiscorr(:,run)/seedvox; % averaging  
                    end
                    
                end
                fX{s}(:,mcfi)=mean(thiscorr,2);

                if writeoutsinglefiles
                    ccmap=msk;
                    ccmap.img=single(ccmap.img);
                    ccmap.fname=[outputfolder,seedfn{s},'_',subIDs{mcfi}{1},'_corr.nii'];
                    ccmap.img(omaskidx)=fX{s}(:,mcfi);    

                    spm_write_vol(ccmap,ccmap.img);
                end
            end
            
            
        otherwise
            for run=1:howmanyruns
                load([dfoldvol,subIDs{mcfi}{run+1}])
                gmtc=single(gmtc);
                
                for s=1:numseed
                    stc(s,:)=mean(gmtc(sweightidx{s},:).*sweightidxmx{s});
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
                save([outputfolder,addp,'corrMx_',subIDs{mcfi}{1},'.mat'],'X','-v7.3');
            end
    end
end
ea_dispercent(1,'end');

switch cmd
    case {'seed','seedvox_ram','seedvox_noram'}
        for s=1:numseed
            % export mean            
            M=nanmean(fX{s}');
            mmap=msk;
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
            
            % lh surf
            lM=nanmean(lh.fX{s}');
            lmmap=lmsk;
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
            rmmap=rmsk;
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
            
            
            % fisher-transform:
            fX{s}=atanh(fX{s});
            lh.fX{s}=atanh(lh.fX{s});
            rh.fX{s}=atanh(rh.fX{s});
            
            % export fz-mean
            
            M=nanmean(fX{s}');
            mmap=msk;
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
            
            % lh surf
            lM=nanmean(lh.fX{s}');
            lmmap=lmsk;
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
            rmmap=rmsk;
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
            
            
            % export T
            
            [~,~,~,tstat]=ttest(fX{s}');
            tmap=msk;
            tmap.img(:)=0;
            tmap.dt=[16,0];
            tmap.img=single(tmap.img);
            
                tmap.img(omaskidx)=tstat.tstat;

            tmap.fname=[outputfolder,seedfn{s},'_func_',cmd,'_T.nii'];
            spm_write_vol(tmap,tmap.img);
            if usegzip
                gzip(tmap.fname);
                delete(mmap.fname);
            end
            
            
            
            
            
            
            % lh surf
            [~,~,~,ltstat]=ttest(lh.fX{s}');
            lmmap=lmsk;
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
            rmmap=rmsk;
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
        
    otherwise
        
        % export mean
        M=nanmean(fX');
        X=zeros(numseed);
        X(logical(triu(ones(numseed),1)))=M;
        X=X+X';
        X(logical(eye(length(X))))=1;
        save([outputfolder,cmd,'_corrMx_AvgR.mat'],'X','-v7.3');
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


toc

function sfile=getrois(sfile)

fID=fopen(sfile);
sfile=textscan(fID,'%s');
sfile=sfile{1};
fclose(fID);







