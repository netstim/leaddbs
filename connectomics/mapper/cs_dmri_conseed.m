function cs_dmri_conseed(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,space)
%# ea_load_nii




[sfile,roilist]=ea_handleseeds(sfile);


if isdeployed
    cbase=dfold;
else
    cbase=ea_getconnectomebase;
end

switch cmd
    case 'seed'
        for s=1:length(sfile)
            map=ea_load_nii([cbase,'spacedefinitions',filesep,space]);
            cfile=[dfold,'dMRI',filesep,cname];
            
            if exist([cfile,filesep,'data.mat'],'file') % regular mat file
                if ~exist('fibers','var')
                    load([cfile,filesep,'data.mat'],'fibers');
                end
                
                redotree=0;
            elseif exist([cfile,filesep,'data.fib.gz'],'file') % regular .fib.gz file
                
                ftr=track_seed_gqi([cfile,filesep,'data.fib.gz'],sfile{s});
                fibers=ftr.fibers;
                redotree=1;
                
            else % connectome type not supported
                ea_error('Connectome file vanished or not supported!');
            end
            
            
            mapsz=size(map.img);
            
            seedfiles=sfile;
            if ~exist('tree','var') || redotree % only compute for first seed.
                tree=KDTreeSearcher(fibers(:,1:3));
            end
            map.img(:)=0;
            Vseed=ea_load_nii(seedfiles{s});
            
            maxdist=abs(Vseed.mat(1))/2;
            
            Vseed.img(isnan(Vseed.img))=0;
            
            ixs=find(Vseed.img);
            % subtract nan values from these
            
            ixvals=Vseed.img(ixs);
            
            [xx,yy,zz]=ind2sub(size(Vseed.img),ixs);
            XYZvx=[xx,yy,zz,ones(length(xx),1)]';
            clear ixs
            XYZmm=Vseed.mat*XYZvx;
            XYZmm=XYZmm(1:3,:)';
            clear Vseed
            
            ids=rangesearch(tree,XYZmm,maxdist,'distance','chebychev');
            % select fibers for each ix
            ea_dispercent(0,'Iterating voxels');
            ixdim=length(ixvals);
            
            for ix=1:ixdim
                % assign fibers on map with this weighted value.
                fibnos=unique(fibers(ids{ix},4));
                
                allfibcs=fibers(ismember(fibers(:,4),fibnos),1:3);
                allfibcs=round(map.mat\[allfibcs,ones(size(allfibcs,1),1)]');
                allfibcs(:,logical(sum(allfibcs<1,1)))=[];
                topaint=sub2ind(mapsz,allfibcs(1,:),allfibcs(2,:),allfibcs(3,:));
                map.img(topaint)=map.img(topaint)+ixvals(ix);
                ea_dispercent(ix/ixdim);
            end
            ea_dispercent(1,'end');
            
            [pth,fn]=fileparts(seedfiles{s});
            
            map.fname=fullfile(outputfolder,[fn,'_struc_',cmd,'.nii']);
            map.dt=[16,0];
            spm_write_vol(map,map.img);
            
        end
    otherwise
        warning('Structural connectivity only supported for seeds');
end




function ftr=track_seed_gqi(cfile,seedfile)


basedir = [ea_getearoot, 'ext_libs',filesep,'dsi_studio',filesep];
if ismac
    dsistudio = [basedir,'mac',filesep, 'dsi_studio.app',filesep,'Contents',filesep,'MacOS',filesep,'dsi_studio'];
elseif isunix
    ea_libs_helper([basedir, 'linux']);
    dsistudio = [basedir, 'linux',filesep,'dsi_studio'];
elseif ispc
    dsistudio = [basedir, 'win',filesep,'dsi_studio.exe'];
end


[pth,fn]=fileparts(seedfile);

cmd=[dsistudio,' --action=trk --source=',ea_path_helper(cfile),...
    ' --method=0',...
    ' --seed=',ea_path_helper(seedfile),...
    ' --seed_count=10000',...
    ' --output=',ea_path_helper([pth,filesep,'temp.mat'])];


err=ea_submitcmd(cmd);
if err
    ea_error(['Fibertracking with dsi_studio failed (error code=',num2str(err),').']);
end

% now store tract in lead-dbs format
ea_dispercent(0,'Converting fibers');
fibinfo=load([pth,filesep,'temp.mat']);
fibers=fibinfo.tracts;
idx=fibinfo.length';
clear fibinfo
fibers=fibers';


fibers(:,1)=78.0-fibers(:,1);
fibers(:,2)=76.0-fibers(:,2);
fibers(:,3)=-50.0+fibers(:,3);

clear length
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'
    ea_dispercent(cnt/length(idx));
    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
ea_dispercent(1,'end');



fibers=[fibers,idxv];

ftr.fourindex=1;
ftr.ea_fibformat='1.0';
ftr.fibers=fibers;
ftr.idx=idx;
delete([pth,filesep,'temp.mat']);

