function varargout=ea_normalize_maget_segment(options)

if ischar(options) % return name of method.
    varargout{1}='MAGeT Brain-like Segmentation/Normalization DISTAL atlas (Chakravarty 2013, Ewert 2016)';
    varargout{2}=1; % dummy output
    varargout{3}=1; % hassettings.
    varargout{4}=1; % is multispectral
    return
end

reforce=0;
atlastouse=options.prefs.machine.normsettings.maget_atlasset; % for now, only the distal atlas is supported!

peerfolders=ea_getmagetpeers(options);

%% step 0: check if all subjects have been processed with an ANTs-based normalization function
ea_magetcheck_norm_peers(options,peerfolders)

subdirec=[options.root,options.patientname,filesep];

%% step 1, warp DISTAL back to each peer brain
earoot=ea_getearoot;
atlasbase=[ea_space(options,'atlases'),atlastouse,filesep];
for peer=1:length(peerfolders)

    peerdirec=[peerfolders{peer},filesep];
    poptions=options;
    [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
    poptions.root=[poptions.root,filesep];
    poptions=ea_assignpretra(poptions);

    %             mkdir([peerdirec,'MAGeT']);
    %             mkdir([peerdirec,'MAGeT',filesep,'atlases']);
    %             mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse]);
    %             mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'lh']);
    %             mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'rh']);
    %             mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'midline']);
    %             mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'mixed']);
    %             mkdir([peerdirec,'MAGeT']);
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
                froms{cnt}=ea_niigz([atlasbase,'lh',filesep,atlases.names{atlas}]);
                sub_tos{cnt}=ea_niigz([satlasbase,'lh',filesep,atlases.names{atlas}]);

                cnt=cnt+1;
            case 2 % RH
                froms{cnt}=ea_niigz([atlasbase,'rh',filesep,atlases.names{atlas}]);
                sub_tos{cnt}=ea_niigz([satlasbase,'rh',filesep,atlases.names{atlas}]);

                cnt=cnt+1;
            case 3 % both RH / LH present
                froms{cnt}=ea_niigz([atlasbase,'lh',filesep,atlases.names{atlas}]);
                sub_tos{cnt}=ea_niigz([satlasbase,'lh',filesep,atlases.names{atlas}]);

                froms{cnt+1}=ea_niigz([atlasbase,'rh',filesep,atlases.names{atlas}]);
                sub_tos{cnt+1}=ea_niigz([satlasbase,'rh',filesep,atlases.names{atlas}]);

                cnt=cnt+2;
            case 4 % Mixed
                froms{cnt}=ea_niigz([atlasbase,'mixed',filesep,atlases.names{atlas}]);
                sub_tos{cnt}=ea_niigz([satlasbase,'mixed',filesep,atlases.names{atlas}]);

                cnt=cnt+1;
            case 5 % Midline
                froms{cnt}=ea_niigz([atlasbase,'midline',filesep,atlases.names{atlas}]);
                sub_tos{cnt}=ea_niigz([satlasbase,'midline',filesep,atlases.names{atlas}]);
                cnt=cnt+1;
        end

    end


    %ea_ants_apply_transforms(poptions,froms,tos,1,[peerdirec,poptions.prefs.prenii_unnormalized]);


    %% step 2, generate warps from MNI via peers to the selected patient brain

    if ~exist([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'2sub.nii.gz'],'file') || reforce
        [~,peerpresentfiles]=ea_assignpretra(poptions);
        [~,subpresentfiles]=ea_assignpretra(options);

        [~,presentinboth]=ismember(subpresentfiles,peerpresentfiles);
        subpresentfiles=subpresentfiles(logical(presentinboth));
        % check the other way:
        [~,presentinboth]=ismember(peerpresentfiles,subpresentfiles);
        peerpresentfiles=peerpresentfiles(logical(presentinboth));

        if ~isequal(subpresentfiles,peerpresentfiles) % then I did something wrong.

        else
            presentfiles=subpresentfiles;
            clear subpresentfiles peerpresentfiles
        end
        if isempty(presentfiles)
            ea_error(['Please supply a peer/subject with valid anatomy files: ',poptions.patientname,'/',options.patientname,'.']);
        end
        clear sptos spfroms weights %clear all variables, if not, the last image file will be transferred onto the next peer
        for anatfi=1:length(presentfiles)
            spfroms{anatfi}=[peerdirec,presentfiles{anatfi}];
            sptos{anatfi}=[subdirec,presentfiles{anatfi}];
        end
        weights=repmat(1.5,length(presentfiles),1);

        % add FA if present ? add to beginning since will use last entry to
        % converge
        if exist([subdirec,options.prefs.fa2anat],'file') && exist([peerdirec,options.prefs.fa2anat],'file')
            spfroms=[{[peerdirec,options.prefs.fa2anat]},spfroms];
            sptos=[{[subdirec,options.prefs.fa2anat]},sptos];

            weights=[0.5;weights];

        end

        defoutput=[subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname];
        if ~exist([subdirec,'MAGeT',filesep,'warpreceives',filesep],'file')
            mkdir([subdirec,'MAGeT',filesep,'warpreceives',filesep]);
        end



        ea_ants_nonlinear(sptos,spfroms,[subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'.nii'],weights,options);
        delete(ea_niigz([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'.nii'])); % we only need the warp
        delete([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'InverseComposite.h5']); % we dont need the inverse warp

        % Now export composite transform from MNI -> Peer -> Subject

        antsdir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
        if ispc
            applyTransforms = ea_path_helper([antsdir, 'antsApplyTransforms.exe']);
        else
            applyTransforms = [antsdir, 'antsApplyTransforms.', computer('arch')];
        end

        prenii=ea_niigz([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
        icmd=[applyTransforms,' -r ',ea_path_helper(prenii),...
              ' -t ',ea_path_helper([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'Composite.h5']),...
              ' -t ',ea_path_helper([peerfolders{peer},filesep,'glanatInverseComposite.h5']),...
              ' -o [',ea_path_helper([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'2sub.nii.gz']),',1]'];
        if ~ispc
            system(['bash -c "', icmd, '"']);
        else
            system(icmd);
        end
        % can cleanup the Peer -> patient transform already.
        delete([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'Composite.h5']); % we dont need the inverse warp
    else
        [~,peerpresentfiles]=ea_assignpretra(poptions);
        [~,subpresentfiles]=ea_assignpretra(options);

        [~,presentinboth]=ismember(subpresentfiles,peerpresentfiles);
        subpresentfiles=subpresentfiles(logical(presentinboth));
        % check the other way:
        [~,presentinboth]=ismember(peerpresentfiles,subpresentfiles);
        presentfiles=peerpresentfiles(logical(presentinboth));

    end


    %% step 3: warp all atlas nuclei from MNI via peer to sub

    if ~exist(sub_tos{end},'file') % assume the work has been done already
        transformfile=[subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'2sub.nii.gz'];
        ea_ants_apply_transforms(poptions,froms,sub_tos,0,[subdirec,presentfiles{1}],transformfile);
    end
    % gather all files for majority voting
    warpednuclei{peer}=sub_tos;
end



%% step 4: Perform majority voting on final atlas
if ~exist([subdirec,'anat_atlas.nii.gz'],'file')
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
        % X(X<0.5)=0; % binarize each image (initial DISTAL atlas is binary, but
        % when warping nonlinearly it can be that two former voxels which intensity value 1 and 0 form only
        % one target voxel with intensitiy value 0.5, (or 0.333 depending on the
        % number of voxels), so this first binarizes the files on the subject
        % level and makes them more conservative and then fuses all images and
        % performs again the majority voting on the overlaid images which also
        % makes the result more conservative)
        % X(X>0)=1;
        X=mean(X,4);

        X=X/max(X(:));
        % X(X<0.25)=0; %defines sensitivitiy of majority voting
        % X(X>0)=1;
        % save atlas
        [pth,atlasname]=fileparts(warpednuclei{peer}{atlas});

        if length(atlasname)>3 && strcmp('.nii',atlasname(end-3:end))
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

    %% step 5: Normalize using multimodal DISTAL warp

    % write out anat_distal
    nii.fname=[subdirec,'anat_atlas.nii'];
    AllX(AllX>1)=1;
    AllX=smooth3(AllX,'gaussian',[5 5 5]);
    nii.img=AllX;
    ea_write_nii(nii);
    gzip(nii.fname);
    delete(nii.fname);
end
if strcmp(options.prefs.dev.profile,'se')
    ; % do siobhanspecific stuff (here do nothing)
else
    ea_normalize_ants(options,1); % do general user specific stuff, here --> performs normalization with
% DISTAL as anchorpoint, commented when segmentation is wanted and not
% normalization, also some results of this normalization are off
end



function ea_writecompositewarp(transforms)
basedir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

cmd = applyTransforms;
refim=[ea_space(options),options.primarytemplate,'.nii'];
% add transforms:
for t=1:length(transforms)
    [pth1,fn1,ext1]=fileparts(transforms{t}{1});
    [pth2,fn2,ext2]=fileparts(transforms{t}{2});
    tr=[' -r ',refim,...
        ' -t [',ea_path_helper([pth1,filesep,fn1,ext1]),',0]',...
        ' -t [',ea_path_helper([pth2,filesep,fn2,ext2]),',0]'];
    cmd = [cmd,tr];
end

% add output:
cmd=[cmd, ' -o [compositeDisplacementField.h5,1]'];

system(cmd)
