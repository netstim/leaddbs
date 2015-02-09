function [DTI_CM] = ea_createCM_dti(options)

maxdist=2;
minlen=3;

if ~exist([options.root,options.patientname,filesep,'w',options.ftCM.atlas,'.nii'],'file')
    %% warp atlas into pre_tra-space:
    
    matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'iy_pre_tra.nii']};
    matlabbatch{1}.spm.util.defs.ofname = '';
    matlabbatch{1}.spm.util.defs.fnames = {[options.earoot,'templates',filesep,'labeling',filesep,options.ftCM.atlas,'.nii,1']};
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.patientname]};
    matlabbatch{1}.spm.util.defs.interp = 1;
    cfg_util('run',{matlabbatch});
    clear matlabbatch
end
if ~exist([options.root,options.patientname,filesep,'rw',options.ftCM.atlas,'.nii'],'file')
    %% coreg atlas into b0-space:
    copyfile([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[options.root,options.patientname,filesep,options.prefs.b0,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[options.root,options.patientname,filesep,'w',options.ftCM.atlas,'.nii,1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    cfg_util('run',{matlabbatch});
    clear matlabbatch
    
    delete([options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    delete([options.root,options.patientname,filesep,'rc',options.prefs.prenii_unnormalized]);
end


%% create voxelmask

Vatl=spm_vol([options.root,options.patientname,filesep,'rw',options.ftCM.atlas,'.nii,1']);
Xatl=spm_read_vols(Vatl);

nonzeros=find(Xatl(:));
vv=Xatl(nonzeros);
[xx,yy,zz]=ind2sub(size(Xatl),nonzeros);

voxelmask.locs=[xx,yy,zz,ones(size(xx,1),1)]';
voxelmask.locs=[Vatl.mat*voxelmask.locs]';
voxelmask.locs=voxelmask.locs(:,1:3);
voxelmask.vals=vv;




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



DTI_CM=zeros(max(voxelmask.vals));





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
Vox=KDTreeSearcher(voxelmask.locs);

%clear seeds terms endpointsearch

    fprintf(1,'Iterating through fibers:         ');


for fiber=1:(fibercount)
    
    percent=round((fiber/fibercount)*100);
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
    
    thisfiblen=length(fibs{fiber});
    thisfib=fibs{fiber}';
    if thisfiblen>minlen % only include fibers >minimum length..
        if sum(thisfib(3,:)<9)<1 % or use isempty(find(thisfib(3,:)<9,1)) % check if fiber exits the brain through spinal chord.. cutoff=9.
            
            
            
            %[sIDX,sD]=knnsearch(Vox,[fibs{fiber}(1,:);fibs{fiber}(end,:)],'IncludeTies',true);
            [sIDX,sD]=rangesearch(Vox,[seeds(fiber,:);terms(fiber,:)],maxdist);

            %seedmin=min(sD(1));
            %targmin=min(sD(4));
            
            
            %if mean(sD{1})<term_mindist && mean(sD{2})<term_mindist % only then make connection.

            if ~isempty(sIDX{1}) && ~isempty(sIDX{2})
            for starts=1:length(sIDX{1})
                for ends=1:length(sIDX{2})
                    weight=1./exp(((sD{1}(starts)+(sD{2}(ends)))/2));
                    
                    
                    DTI_CM(voxelmask.vals(sIDX{1}(starts)),voxelmask.vals(sIDX{2}(ends)))    =  ...
                        DTI_CM(voxelmask.vals(sIDX{1}(starts)),voxelmask.vals(sIDX{2}(ends)))    +   weight; %*thisfiblen; %% +1 is to avoid dividing by zero..
                    DTI_CM(voxelmask.vals(sIDX{2}(ends)),voxelmask.vals(sIDX{1}(starts)))    =  ...
                        DTI_CM(voxelmask.vals(sIDX{2}(ends)),voxelmask.vals(sIDX{1}(starts)))    +   weight; %*thisfiblen;   % symmetrize Matrix.
                    
                end
            end
            conns=conns+1;
            end
        end
    end
    %end
end

     fprintf('\n')


disp(['In total used ',num2str(conns),'/',num2str(fiber),' fibers to connect ',num2str(length(DTI_CM)),' regions.']);

disp('Done.');
fprintf('\n');


save([options.root,options.patientname,filesep,'DTI_CM.mat'],'DTI_CM','-v7.3');

cm=figure('color','w','NumberTitle','off','Name','DTI Connectivity Matrix');
imagesc(DTI_CM);
colorbar
saveas(cm,[options.root,options.patientname,filesep,'DTI_CM.png']);
