function [DTI_CM] = ea_createCM_dti(options)
% This function creates a structural (dMRI-based) Connectivity matrix.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

minlen=options.prefs.lc.struc.minlen;
directory=[options.root,options.patientname,filesep];

%% get node definition of current parcellation scheme
ea_warp_parcellation(options.prefs.b0,'b0',options);
Vatl=spm_vol([directory,'templates',filesep,'labeling',filesep,'rb0w',options.lc.general.parcellation,'.nii,1']);

%% get fiber definition
disp('Loading FTR-File.');
fs=load([[options.root,options.patientname,filesep],options.prefs.FTR_unnormalized]);

if ~isfield(fs,'curveSegCell') % Freiburg Format
    fn = fieldnames(fs);
    eval(sprintf('curveSegCell = fs.%s;',fn{1}));
    if size(curveSegCell,1)>size(curveSegCell,2)
        curveSegCell=curveSegCell';
    end
    convertfromfreiburg=0;
    fibs=curveSegCell;
else
    convertfromfreiburg=1;
    fibs=fs.curveSegCell;
end

clear curveSegCell fs

%ea_dispercent(0,'Calculating seeds and terminals...');
disp('Calculating seeds and terminals...');

fibercount=length(fibs);
seeds=zeros(fibercount,3);
terms=zeros(fibercount,3);
for fiber=1:(fibercount)
 %   ea_dispercent(fiber/fibercount);
    seeds(fiber,:)=fibs{fiber}(1,:);
    terms(fiber,:)=fibs{fiber}(end,:);
end
%ea_dispercent(100,'end');


if convertfromfreiburg % already in voxel notation
   Xatl=spm_read_vols(Vatl);
   ysize=size(Xatl,2)+1;
   seeds=[ysize-seeds(:,2),seeds(:,1),seeds(:,3),ones(size(seeds,1),1)]; % yflip, switch x and y (reversing freiburg notation)
   terms=[ysize-terms(:,2),terms(:,1),terms(:,3),ones(size(terms,1),1)];
   seeds=seeds(:,1:3);
   terms=terms(:,1:3);
else
   if any(seeds(:)<-1) || any(terms(:)<-1) % mm-notation, convert to voxel notation
       seeds=[seeds(:,1),seeds(:,2),seeds(:,3),ones(size(seeds,1),1)]';
       terms=[terms(:,1),terms(:,2),terms(:,3),ones(size(terms,1),1)]';
       seeds=(Vatl.mat\seeds)';
       terms=(Vatl.mat\terms)';
       seeds=seeds(:,1:3);
       terms=terms(:,1:3);
   end
end

seedIDX=int16(spm_sample_vol(Vatl,seeds(:,1), seeds(:,2), seeds(:,3),0));
termIDX=int16(spm_sample_vol(Vatl,terms(:,1), terms(:,2), terms(:,3),0));

clear seeds terms

%% create CM
display('Initializing structural CM.');

aID = fopen([options.earoot,'templates',filesep,'labeling',filesep,options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
d=length(atlas_lgnd{1}); % how many ROI.
DTI_CM=zeros(d);

ea_dispercent(0,['Iterating through ',num2str(length(fibs)),' fibers']);
conns=0;
endpts=[seedIDX,termIDX];
conIDX=endpts>0;
conIDX=sum(conIDX,2);
conIDX=find(conIDX==2); % only take fibers into account that dont start/end in zero-regions
for fiber=conIDX'
    percent=fiber/fibercount;
    ea_dispercent(percent);
    if length(fibs{fiber})>minlen % only include fibers >minimum length
        %if sum(thisfib(3,:)<-55)<1 % check if fiber exits the brain through spinal chord.. choose a cutoff(i.e.=-55mm).
        DTI_CM(seedIDX(fiber),termIDX(fiber))    =  ...
            DTI_CM(seedIDX(fiber),termIDX(fiber))    +   1;
        DTI_CM(termIDX(fiber),seedIDX(fiber))    =  ...
            DTI_CM(seedIDX(fiber),termIDX(fiber));  % symmetrize Matrix.
        conns=conns+1; % connection count
        %end
    end
end

ea_dispercent(100,'end')


disp(['In total used ',num2str(conns),'/',num2str(fiber),' fibers to connect ',num2str(length(DTI_CM)),' regions.']);

disp('Done.');
