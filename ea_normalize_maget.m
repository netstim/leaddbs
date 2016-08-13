function varargout=ea_normalize_maget(options)

if ischar(options) % return name of method.
    varargout{1}='MAGeT Brain-like Segmentation/Normalization DISTAL atlas (Chakravarty 2013, Ewert 2016)';
    varargout{2}={'SPM8','SPM12'};
    return
end

reforce=0;
atlastouse='DISTAL_manual'; % for now, only the distal atlas is supported!

ptage=ea_getpatientage([options.root,options.patientname,filesep]);

peerfolders=ea_getIXI_IDs(3,ptage);

%% step 0: check if all subjects have been processed with an ANTs-based normalization function
for peer=1:length(peerfolders)
   if ~ismember(ea_whichnormmethod([peerfolders{peer},filesep]),ea_getantsnormfuns)
      ea_error('Please make sure that all peers selected have been normalized using ANTs.')
   end
end

subdirec=[options.root,options.patientname,filesep];
if ~ismember(ea_whichnormmethod(subdirec),ea_getantsnormfuns)
   ea_normalize_ants_multimodal(options,'coregonly'); 
end

%% step 1, setup DISTAL warps back to sub via each peer brain
earoot=ea_getearoot;
atlasbase=[earoot,'atlases',filesep,atlastouse,filesep];
for peer=1:length(peerfolders)

    peerdirec=[peerfolders{peer},filesep];
    poptions=options;
    [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
    poptions.root=[poptions.root,filesep];
    poptions=ea_assignpretra(poptions);

            mkdir([subdirec,'MAGeT',filesep,'atlasreceives']);
            mkdir([subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname]);
            mkdir([subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse]);
            mkdir([subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'lh']);
            mkdir([subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'rh']);
            mkdir([subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'midline']);
            mkdir([subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'mixed']);
        
        satlasbase=[subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep];

        load([atlasbase,'atlas_index.mat']);
        cnt=1;
        for atlas=1:length(atlases.names)
            % warp atlas to peer

            switch atlases.types(atlas)
                case 1 % LH only
                    froms{cnt}=[atlasbase,'lh',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'lh',filesep,atlases.names{atlas}];
                    
                    cnt=cnt+1;
                case 2 % RH
                    froms{cnt}=[atlasbase,'rh',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'rh',filesep,atlases.names{atlas}];

                    cnt=cnt+1;
                case 3 % both RH / LH present
                    froms{cnt}=[atlasbase,'lh',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'lh',filesep,atlases.names{atlas}];

                    froms{cnt+1}=[atlasbase,'rh',filesep,atlases.names{atlas}];
                    sub_tos{cnt+1}=[satlasbase,'rh',filesep,atlases.names{atlas}];

                    cnt=cnt+2;
                case 4 % Mixed
                    froms{cnt}=[atlasbase,'mixed',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'mixed',filesep,atlases.names{atlas}];

                    cnt=cnt+1;
                case 5 % Midline
                    froms{cnt}=[atlasbase,'midline',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'midline',filesep,atlases.names{atlas}];

                    cnt=cnt+1;
            end
            
        end
        
        % check for .gz files:
        for fi=1:length(froms)
            if ~exist(froms{fi},'file') && exist([froms{fi},'.gz'],'file')
                froms{fi}=[froms{fi},'.gz'];
                sub_tos{fi}=[sub_tos{fi},'.gz'];
            else
                ea_error(['Atlas file not found or dual .nii.gz/.nii version present: ',froms{fi}]);
            end
        end
            
    
    
    %% step 2, generate warps from peers to the selected patient brain

    
    if ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2mni.nii.gz'],'file') || reforce
        [~,peerpresentfiles]=ea_assignpretra(poptions);
        [~,subpresentfiles]=ea_assignpretra(options);
        presentinboth=ismember(subpresentfiles,peerpresentfiles);
        peerpresentfiles=peerpresentfiles(presentinboth);
        if ~isequal(subpresentfiles,peerpresentfiles) % then I made something wrong.
            keyboard
        else
           presentfiles=subpresentfiles;
           clear subpresentfiles peerpresentfiles
        end
        for anatfi=1:length(presentfiles)
            spfroms{anatfi}=[peerdirec,presentfiles{anatfi}];
            sptos{anatfi}=[subdirec,presentfiles{anatfi}];
            metrics{anatfi}='MI';
        end
        weights=repmat(1.5,length(presentfiles),1);
       
        % add FA if present ? add to beginning since will use last entry to
        % converge
        if exist([subdirec,options.prefs.fa2anat],'file') && exist([peerdirec,options.prefs.fa2anat],'file')
            spfroms=[{[peerdirec,options.prefs.fa2anat]},spfroms];
            sptos=[{[subdirec,options.prefs.fa2anat]},sptos];
            
           weights=[0.5;weights];
           metrics=[{'MI'},metrics];
           
        end
        
        defoutput=[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname];
        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep],'file')
            mkdir([subdirec,'MAGeT',filesep,'warps',filesep]);
        end
        
        
        ea_ants_nonlinear(sptos,spfroms,[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'.nii'],weights,metrics,options);
        delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'.nii']); % we only need the warp
        %delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'InverseComposite.h5']); % we dont need the inverse warp
  
        % Now export composite transform from MNI -> Peer -> Subject
    
        if ispc
            sufx='.exe';
        else
            sufx=computer('arch');
        end
        
        antsApply=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep,'antsApplyTransforms.',sufx];
        
        template=[ea_getearoot,'templates',filesep,'mni_hires.nii'];
        prenii=[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized];
        cmd=[antsApply,' -r ',template,' -t ',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'InverseComposite.h5'],' -t ',[peerfolders{peer},filesep,'glanatInverseComposite.h5'],' -o [',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2mni.nii.gz',',1]']];
        icmd=[antsApply,' -r ',template,' -t ',[peerfolders{peer},filesep,'glanatInverseComposite.h5'],' -t ',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5'],' -o [',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2sub.nii.gz',',1]']];
            if ~ispc
                system(['bash -c "', cmd, '"']);
                system(['bash -c "', icmd, '"']);
            else
                system(cmd);
                system(icmd);
            end
            
            % delete intermediary files
            delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'InverseComposite.h5']);
            delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5']);   
    end
    
    
    %% step 3: warp all atlas nuclei from peers to sub
    
    %if ~exist([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse],'file') % assume the work has been done already
        % caution: old tos is new from here! this is correct.
        transformfile=[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2sub.nii.gz'];
        ea_ants_applytransforms(poptions,froms,sub_tos,0,[subdirec,poptions.prefs.prenii_unnormalized],transformfile);
        
    %end

    % gather all files for majority voting
    warpednuclei{peer}=sub_tos;
end


% aggregate all warps together and average them

if ispc
    sufx='.exe';
else
    sufx=computer('arch');
end
antsAve=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep,'antsAverageImages.',sufx];
warpbase=[options.root,options.patientname,filesep,'MAGeT',filesep,'warps'];
cmd=[antsAve,' 3 ',warpbase,'ave2mni.nii.gz',' 0 ',warpbase,'*2mni.nii.gz'];
icmd=[antsAve,' 3 ',warpbase,'ave2sub.nii.gz',' 0 ',warpbase,'*2sub.nii.gz'];
if ~ispc
    system(['bash -c "', cmd, '"']);
    system(['bash -c "', icmd, '"']);
else
    system(cmd);
    system(icmd);
end


keyboard

%% step 4: Perform majority voting on final atlas
mkdir([subdirec,'atlases']);
mkdir([subdirec,'atlases',filesep,'native']);
mkdir([subdirec,'atlases',filesep,'native',filesep,atlastouse]);
mkdir([subdirec,'atlases',filesep,'native',filesep,atlastouse,filesep,'lh']);
mkdir([subdirec,'atlases',filesep,'native',filesep,atlastouse,filesep,'rh']);
mkdir([subdirec,'atlases',filesep,'native',filesep,atlastouse,filesep,'mixed']);
mkdir([subdirec,'atlases',filesep,'native',filesep,atlastouse,filesep,'midline']);
for atlas=1:length(warpednuclei{peer})
    for peer=1:length(peerfolders) % seek to average across peers
        nii=ea_load_nii(warpednuclei{peer}{atlas}); % add gz support
        if peer==1
            X=zeros([size(nii.img),length(warpednuclei{peer})]);
        end
        X(:,:,:,peer)=nii.img;
    end
    X(X<0.5)=0; % binarize each image
    X(X>0)=1;
    X=mean(X,4);
    
    X=X/max(X(:));
    X(X<0.5)=0;
    X(X>0)=1;
    % save atlas
    [pth,atlasname]=fileparts(warpednuclei{peer}{atlas});
    if strcmp('.nii',atlasname(end-3:end))
        atlasname=atlasname(1:end-4);
    end
    [~,base]=fileparts(pth);
    nii.fname=[subdirec,'atlases',filesep,'native',filesep,atlastouse,filesep,base,filesep,atlasname,'.nii'];
    nii.img=X;
    if ~exist('AllX','var')
        AllX=X;
    else
        AllX=AllX+X;
    end
    
    clear X
    ea_write_nii(nii);
    ea_crop_nii(nii.fname);
    gzip(nii.fname);
    delete(nii.fname);
end


% -> Now segmentation is done. Add normalization.
switch normmethod
    case 'atlasbased'
%% step 5: Normalize using multimodal DISTAL warp

% write out anat_distal
nii.fname=[subdirec,'anat_atlas.nii'];
AllX(AllX>1)=1;
AllX=smooth3(AllX,'gaussian',[5 5 5]);
nii.img=AllX;
ea_write_nii(nii);
gzip(nii.fname);
delete(nii.fname);

ea_normalize_ants_multimodal(options,0,1); 

    case 'wholebrain'
        
        
end

function ea_writecompositewarp(transforms)
basedir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
if ispc
    applyTransforms = [basedir, 'antsApplyTransforms.exe'];
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end


cmd=[applyTransforms];
refim=[options.earoot,'templates',filesep,'mni_hires.nii'];
% add transforms:
for t=1:length(transforms)
    [pth1,fn1,ext1]=fileparts(transforms{t}{1});
    [pth2,fn2,ext2]=fileparts(transforms{t}{2});
        tr=[' -r ',refim,...
        ' -t [',ea_path_helper([pth1,filesep,fn1,ext1]),',0]',...
        ' -t [',ea_path_helper([pth2,filesep,fn2,ext2]),',0]'];
cmd=[cmd,tr];
end

% add output:
cmd=[cmd,' -o [compositeDisplacementField.h5,1]'];

system(cmd)

