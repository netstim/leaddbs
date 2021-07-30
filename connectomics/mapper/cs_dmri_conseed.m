function cs_dmri_conseed(connBaseFolder, options)

connName = options.lcm.struc.connectome;
sfile = ea_handleseeds(options.lcm.seeds');
cmd = ea_lcm_resolvecmd(options.lcm.cmd);

switch options.lcm.struc.espace
    case 1
        space = '222.nii';
    case 2
        space = '111.nii';
    case 3
        space = '555.nii';
end

if strcmp(cmd, 'seed') % Calculate seed connectivity, output folder will be derived automatically
	outputfolder=[];
else
    outputfolder = options.lcm.odir;
    if isempty(outputfolder)
        warning('off', 'backtrace');
        warning('Custom output folder not specified! Will save result to current folder.');
        warning('on', 'backtrace');
        outputfolder = ea_getoutputfolder({[pwd, filesep]},connName);
    elseif ~strcmp(outputfolder(end),filesep)
        outputfolder = [outputfolder,filesep];
    end
end

disp(['Command: ',cmd]);
switch cmd
    case 'seed'             
        cs_dmri_conseed_map(connBaseFolder,connName,sfile,cmd,space,options)
    case {'matrix', 'pmatrix'}
        for s=1:length(sfile)
            if strcmp(connBaseFolder, 'Patient''s fiber tracts')
                if strcmp(connName, options.prefs.FTR_normalized) % patient specific fibertracts
                    cfile=[options.uivatdirs{s},filesep,'connectomes',filesep,'dMRI',filesep,'wFTR.mat'];
                    [fibers,fidx]=ea_loadfibertracts(cfile);
                    redotree=1;
                else % connectome type not supported
                    ea_error(['Connectome file (',options.prefs.FTR_normalized,') vanished or not supported!']);
                end
            else
                cfile=[connBaseFolder,'dMRI',filesep,connName];

                if exist([cfile,filesep,'data.mat'],'file') % regular mat file
                    if ~exist('fibers','var')
                        [fibers,fidx]=ea_loadfibertracts([cfile,filesep,'data.mat']);
                        if ~exist('fibers','var')
                            ea_error('Structural connectome file supplied in wrong format.');
                        end
                    end
                    redotree=0;
                elseif exist([cfile,filesep,'data.fib.gz'],'file') % regular .fib.gz file
                    ftr=track_seed_gqi([cfile,filesep,'data.fib.gz'],sfile{s});
                    fibers=ftr.fibers;
                    redotree=1;
                else % connectome type not supported
                    ea_error('Connectome file vanished or not supported!');
                end
            end

            Vseed=ea_load_nii(sfile{s});

            maxdist=mean(abs(Vseed.voxsize))/2;

            Vseed.img(isnan(Vseed.img))=0;

            ixs=find(Vseed.img);
            % subtract nan values from these

            ixvals=Vseed.img(ixs);

            % now have all seeds and connectome - for each seed find fibers
            % connected to it:

            [xx,yy,zz]=ind2sub(size(Vseed.img),ixs);
            XYZvx=[xx,yy,zz,ones(length(xx),1)]';
            XYZmm=Vseed.mat*XYZvx;
            XYZmm=XYZmm(1:3,:)';
            if ~exist('tree','var') || redotree % only compute for first seed.
                tree=KDTreeSearcher(fibers(:,1:3));
            end
            ids=rangesearch(tree,XYZmm,maxdist,'distance','chebychev');

            % select fibers for each ix
            ea_dispercent(0, ['Iterating voxels (', num2str(s), '/', num2str(length(sfile)), ')']);
            fiberstrength{s}=zeros(size(fidx,1),1); % in this var we will store a mean value for each fiber (not fiber segment) traversing through seed
            fiberstrengthn{s}=zeros(size(fidx,1),1); % counting variable to average strengths
            for ix=find(cellfun(@length,ids))' % iterate through voxels that found something
                % assign fibers on map with this weighted value.
                fibnos=unique(fibers(ids{ix},4)); % these fiber ids go through this particular voxel.
                fiberstrength{s}(fibnos)=fiberstrength{s}(fibnos)+ixvals(ix);
                fiberstrengthn{s}(fibnos)=fiberstrengthn{s}(fibnos)+1;
                % ea_dispercent(ix/ixdim);
            end
            nzz=~(fiberstrength{s}==0);
            fiberstrength{s}(nzz)=fiberstrength{s}(nzz)./fiberstrengthn{s}(nzz); % now each fiber has a strength mediated by the seed.
            ea_dispercent(1,'end');
        end

        fiberstrength=cell2mat(fiberstrength);
        mat=zeros(length(sfile));

        for sxx=1:length(sfile)
            for syy=1:length(sfile)
                if sxx>syy
                    switch cmd
                        case 'matrix'
                            mask=logical((fiberstrength(:,sxx)>0).*(fiberstrength(:,syy)>0));
                        case 'pmatrix'
                            oix=1:length(sfile);
                            oix([sxx,syy])=[];
                            mask=logical((fiberstrength(:,sxx)>0).*(fiberstrength(:,syy)>0).*(~(any(fiberstrength(:,oix)>0,2))));
                    end
                    mat(sxx,syy)=sum(fiberstrength(mask,sxx)+fiberstrength(mask,syy))/2;

                elseif sxx==syy % also fill diagonal
                    mat(sxx,syy)=sum(fiberstrength(:,sxx));
                end
            end
        end
        mat=mat+mat'; % symmetrize matrix
        mat(logical(eye(length(sfile))))=mat(logical(eye(length(sfile))))/2;
        if isfield(options.lcm, 'parcSeedName') && ~isempty(options.lcm.parcSeedName)
            fn = options.lcm.parcSeedName;
        else
            [~,fn]=fileparts(sfile{1});
        end
        save(fullfile(outputfolder,[fn,'_struc_',cmd,'.mat']),'mat');
    otherwise
        warning('Structural connectivity only supported for seed / matrix / pmatrix commands.');
end


function ftr=track_seed_gqi(cfile,seedfile)

basedir = [ea_getearoot, 'ext_libs',filesep,'dsi_studio',filesep];
if ispc
    DSISTUDIO = ea_path_helper([basedir, 'dsi_studio.exe']);
else
    DSISTUDIO = [basedir, 'dsi_studio.', computer('arch')];
end

pth=fileparts(seedfile);

cmd=[DSISTUDIO,' --action=trk --source=',ea_path_helper(cfile),...
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
