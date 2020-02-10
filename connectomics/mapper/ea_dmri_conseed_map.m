function ea_dmri_conseed_map(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,space,options)

owasempty=0;
if ~exist('outputfolder','var')
    outputfolder=ea_getoutputfolder(sfile,ocname);
    owasempty=1;
else
    if isempty(outputfolder) % from shell wrapper.
        outputfolder=ea_getoutputfolder(sfile,cname);
        owasempty=1;
    end
    if ~strcmp(outputfolder(end),filesep)
        outputfolder=[outputfolder,filesep];
    end
end

for s=1:length(sfile)
    
    % okay to check to reload connectome again for each seed in case of
    % using patient-specific connectomes:
    map=ea_load_nii([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,space]);
    
    if strcmp(dfold, 'Patient''s fiber tracts')
        if strcmp(cname, options.prefs.FTR_normalized) % patient specific fibertracts
            cfile=[options.uivatdirs{s},filesep,'connectomes',filesep,'dMRI',filesep,'wFTR.mat'];
            [fibers,fidx,voxmm,mat]=ea_loadfibertracts([cfile]);
            redotree=1;
            ctype='mat';
        else % connectome type not supported
            ea_error(['Connectome file (',options.prefs.FTR_normalized,') vanished or not supported!']);
        end
    else
        cfile=[dfold,'dMRI',filesep,cname];
        if exist([cfile,filesep,'data.mat'],'file') % regular mat file
            if ~exist('fibers','var')
                [fibers,fidx,voxmm,mat]=ea_loadfibertracts([cfile,filesep,'data.mat']);
                if ~exist('fibers','var')
                    ea_error('Structural connectome file supplied in wrong format.');
                end
            end
            
            redotree=0;
            ctype='mat';
        elseif exist([cfile,filesep,'data.fib.gz'],'file') % regular .fib.gz file
            ftr=track_seed_gqi([cfile,filesep,'data.fib.gz'],sfile{s});
            fibers=ftr.fibers;
            redotree=0;
            ctype='fibgz';
        else % connectome type not supported
            ea_error('Connectome file vanished or not supported!');
        end
    end
    
    mapsz=size(map.img);
    
    Niter=1; % usually run seed only once
    if evalin('base','exist(''SB_SEED_BOUNCE'',''var'')')
        Niter=10;
    end
    
    for iter=1:Niter
        if iter>1 % bounce iteration
            Vseed=map;
%            Vseed.img(~(Vseed.img==0))=ea_minmax(Vseed.img(~(Vseed.img==0))); % reset scale from 0-1 % need to think about negatives
        else
            Vseed=ea_load_nii(sfile{s});
        end
        map.img(:)=0;

        
        maxdist=mean(abs(Vseed.voxsize))/2;
        
        Vseed.img(isnan(Vseed.img))=0;
        
        ixs=find(Vseed.img);
        % subtract nan values from these
        
        ixvals=Vseed.img(ixs);
        if sum(abs(ixvals-double(logical(ixvals))))<0.0001
            allbinary=1;
        else
            allbinary=0;
        end
        
        if ~allbinary || strcmp(ctype,'mat')
            
            [xx,yy,zz]=ind2sub(size(Vseed.img),ixs);
            XYZvx=[xx,yy,zz,ones(length(xx),1)]';
            clear ixs
            XYZmm=Vseed.mat*XYZvx;
            XYZmm=XYZmm(1:3,:)';
            %clear Vseed
            if ~exist('tree','var') || redotree % only compute for first seed.
                tree=KDTreeSearcher(fibers(:,1:3));
            end
            ids=rangesearch(tree,XYZmm,maxdist,'distance','chebychev');
            % select fibers for each ix
            ea_dispercent(0, ['Iterating seeds (', num2str(s), '/', num2str(length(sfile)), ')']);
            ixdim=length(ixvals);
            fiberstrength=zeros(size(fidx,1),1); % in this var we will store a mean value for each fiber (not fiber segment) traversing through seed
            fiberstrengthn=zeros(size(fidx,1),1); % counting variable to average strengths
            for ix=find(cellfun(@length,ids))'
                % assign fibers on map with this weighted value.
                %if ~isempty(ids{ix})
                fibnos=unique(fibers(ids{ix},4)); % these fiber ids go through this particular voxel.
                fiberstrength(fibnos)=fiberstrength(fibnos)+ixvals(ix);
                fiberstrengthn(fibnos)=fiberstrengthn(fibnos)+1;
                %end
                % ea_dispercent(ix/ixdim);
                
            end
            nzz=~(fiberstrength==0);
            fiberstrength(nzz)=fiberstrength(nzz)./fiberstrengthn(nzz); % now each fiber has a strength mediated by the seed.
            ea_dispercent(1,'end');
            
            ea_dispercent(0, ['Iterating fibers (', num2str(s), '/', num2str(length(sfile)), ')']);
            cfibers=find(fiberstrength);
            
            allfibcs = fibers(ismember(fibers(:,4),cfibers), 1:3);
            allfibcs = round(map.mat\[allfibcs, ones(size(allfibcs,1),1)]');
            todel = logical(sum(allfibcs<1,1));
            allfibcs(:, todel) = [];
            topaint = sub2ind(mapsz, allfibcs(1,:), allfibcs(2,:), allfibcs(3,:));
            
            fibInd = fibers(ismember(fibers(:,4),cfibers), 4)';
            fibInd(todel) = [];
            topaint = splitapply(@(x) {unique(x)}, topaint, findgroups(fibInd));
            
            fibInd = repelem(unique(fibInd,'stable'), cellfun(@length, topaint));
            topaint = cell2mat(topaint);
            
            [uniqueImgInd, ~, ic] = unique(topaint);
            fiberStr = accumarray(ic, fiberstrength(fibInd))';
            
            map.img(uniqueImgInd) = map.img(uniqueImgInd) + fiberStr;
            
            ea_dispercent(1,'end');
        else % if all is binary && using a .fib.gz file (i.e. all fibers go through seed already), can be much quicker.
            allfibcs=fibers(:,1:3);
            allfibcs=round(map.mat\[allfibcs,ones(size(allfibcs,1),1)]');
            topaint=sub2ind(mapsz,allfibcs(1,:),allfibcs(2,:),allfibcs(3,:));
            utopaint=unique(topaint);
            c=countmember(utopaint,topaint);
            map.img(utopaint)=c;
        end
        
        [~,fn]=fileparts(sfile{s});
        if owasempty
            outputfolder=ea_getoutputfolder({sfile{s}},cname);
        end
        
        if evalin('base','exist(''SB_SEED_BOUNCE'',''var'')')
            map.img(~(map.img==0))=ea_normal(map.img(~(map.img==0)));
        end
        
        map.fname=fullfile(outputfolder,[fn,'_struc_',cmd,'.nii']);
        map.dt=[16,0];
        spm_write_vol(map,map.img);
        if evalin('base','exist(''SB_SAVE_ITERS'',''var'')')
            copyfile(fullfile(outputfolder,[fn,'_struc_',cmd,'.nii']),fullfile(outputfolder,[fn,'_struc_',cmd,'_',num2str(iter),'.nii']));
        end
        
        
        if iter>1
           disp(['Similarity to last: ',num2str(corr(map.img(:),Vseed.img(:))),'.']); 
        end
    end
end

function D=ea_minmax(D)

D(:)=D(:)-(min(D(:)));
D(:)=D(:)./max(D(:));
