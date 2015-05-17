function [DTI_CM] = ea_createCM_dti(options)

maxdist=options.lc.struc.maxdist;
minlen=options.prefs.lc.struc.minlen;

directory=[options.root,options.patientname,filesep];

ea_warp_parcellation(options.prefs.b0,'b0',options);


%% create voxelmask

Vatl=spm_vol([directory,'templates',filesep,'labeling',filesep,'rb0w',options.lc.general.parcellation,'.nii,1']);
Xatl=spm_read_vols(Vatl);

nonzeros=find(Xatl(:));
vv=Xatl(nonzeros);
[xx,yy,zz]=ind2sub(size(Xatl),nonzeros);

voxelmask.locsvx=[xx,yy,zz,ones(size(xx,1),1)]';
voxelmask.locsmm=[Vatl.mat*voxelmask.locsvx]';
voxelmask.locsvx=voxelmask.locsvx(:,1:3);
voxelmask.locsmm=voxelmask.locsmm(:,1:3);
voxelmask.vals=round(vv);





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




display('Initializing Structural CM.');

aID = fopen([options.earoot,'templates',filesep,'labeling',filesep,options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
d=length(atlas_lgnd{1}); % how many ROI.
DTI_CM=zeros(d);





disp(['Iterating through ',num2str(length(fibs)),' fibers and ',num2str(length(voxelmask.vals)),' Voxels.']);
disp(['Maximum distance for connections is ',num2str(maxdist),'.']);

fibercount=length(fibs);
conns=0;



seeds=zeros(fibercount,3);
terms=zeros(fibercount,3);

disp('Calculating seeds and terminals...');
for fiber=1:(fibercount)
    
    seeds(fiber,:)=fibs{fiber}(1,:);
    terms(fiber,:)=fibs{fiber}(end,:);
    %endpoints=[seeds;terms];
end





if convertfromfreiburg
   Vb0=spm_vol([options.root,options.patientname,filesep,options.prefs.b0]);
   Xb0=spm_read_vols(Vb0);
   ysize=size(Xb0,2)+1;
   seeds=[ysize-seeds(:,2),seeds(:,1),seeds(:,3),ones(size(seeds,1),1)]'; % yflip, switch x and y (reversing freiburg notation)
   terms=[ysize-terms(:,2),terms(:,1),terms(:,3),ones(size(terms,1),1)]';
   
   seeds=[Vb0.mat*seeds]';
   terms=[Vb0.mat*terms]';
   
   seeds=seeds(:,1:3);
   terms=terms(:,1:3);
else
   if ~any(seeds(:)<-1) && ~any(terms(:)<-1) % this would make (erroneous) voxel-notation highly probable..
      warning('It seems like fibers have been stored in voxel dimensions. Converting to mm notation...'); 
      seeds=[seeds(:,1),seeds(:,2),seeds(:,3),ones(size(seeds,1),1)]';
      terms=[terms(:,1),terms(:,2),terms(:,3),ones(size(terms,1),1)]';
      
      seeds=[Vb0.mat*seeds]';
      terms=[Vb0.mat*terms]';
      
      seeds=seeds(:,1:3);
      terms=terms(:,1:3);
      
   end
    
end


disp('Establishing treesearcher for voxelmask...');
Vox=KDTreeSearcher(voxelmask.locsmm);

%clear seeds terms endpointsearch

    fprintf(1,'Iterating through fibers:         ');


for fiber=1:(fibercount)
    
    percent=round((fiber/fibercount)*100);
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
    
    thisfiblen=length(fibs{fiber});
    %thisfib=fibs{fiber}';
    if thisfiblen>minlen % only include fibers >minimum length..
        %if sum(thisfib(3,:)<-55)<1 % check if fiber exits the brain through spinal chord.. cutoff=9.
            
            
            
            sIDX=rangesearch(Vox,[seeds(fiber,:);terms(fiber,:)],maxdist/2,'distance','chebychev');

            %seedmin=min(sD(1));
            %targmin=min(sD(4));
            
            
            %if mean(sD{1})<term_mindist && mean(sD{2})<term_mindist % only then make connection.

            if ~isempty(sIDX{1}) && ~isempty(sIDX{2})
            for starts=1:length(sIDX{1})
                for ends=1:length(sIDX{2})
                    %weight=1; %./exp(((sD{1}(starts)+(sD{2}(ends)))/2));
                    
                    
                    DTI_CM(voxelmask.vals(sIDX{1}(starts)),voxelmask.vals(sIDX{2}(ends)))    =  ...
                        DTI_CM(voxelmask.vals(sIDX{1}(starts)),voxelmask.vals(sIDX{2}(ends)))    +   1; %*thisfiblen; %% +1 is to avoid dividing by zero..
                    DTI_CM(voxelmask.vals(sIDX{2}(ends)),voxelmask.vals(sIDX{1}(starts)))    =  ...
                        DTI_CM(voxelmask.vals(sIDX{2}(ends)),voxelmask.vals(sIDX{1}(starts)))    +   1; %*thisfiblen;   % symmetrize Matrix.
                    
                end
            end
            conns=conns+1;
            end
        %end
    end
    %end
end

     fprintf('\n')


disp(['In total used ',num2str(conns),'/',num2str(fiber),' fibers to connect ',num2str(length(DTI_CM)),' regions.']);

disp('Done.');
fprintf('\n');
