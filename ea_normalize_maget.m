function varargout=ea_normalize_maget(options)

if ischar(options) % return name of method.
    varargout{1}='MAGeT Brain-like Normalization (Chakravarty 2013)';
    varargout{2}=1; % is compatible
    varargout{3}=1; % hassettings.
    varargout{4}=1; % is multispectral
    return
end

reforce=1;

peerfolders=ea_getmagetpeers(options);


%% step 0: check if all subjects have been processed with an ANTs-based normalization function
ea_magetcheck_norm_peers(options,peerfolders)


subdirec=[options.root,options.patientname,filesep];


%% step 1, setup DISTAL warps back to sub via each peer brain
earoot=ea_getearoot;

for peer=1:length(peerfolders)

    clear spfroms sptos weights
    peerdirec=[peerfolders{peer},filesep];
    poptions=options;
    [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
    poptions.root=[poptions.root,filesep];
    poptions=ea_assignpretra(poptions);


    %% step 1, generate warps from peers to the selected patient brain


    if ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2mni.nii'],'file') || reforce
        [~,peerpresentfiles]=ea_assignpretra(poptions);
        [~,subpresentfiles]=ea_assignpretra(options);
        [~,presentinboth]=ismember(subpresentfiles,peerpresentfiles);
        subpresentfiles=subpresentfiles(logical(presentinboth));
        % check the other way:
        [~,presentinboth]=ismember(peerpresentfiles,subpresentfiles);
        peerpresentfiles=peerpresentfiles(logical(presentinboth));
        if ~isequal(subpresentfiles,peerpresentfiles) % then I did something wrong.
            keyboard
        else
            presentfiles=subpresentfiles;
            clear subpresentfiles peerpresentfiles
        end
        clear sptos spfroms weights
        for anatfi=1:length(presentfiles)
            spfroms{anatfi}=ea_niigz([subdirec,presentfiles{anatfi}]);
            sptos{anatfi}=ea_niigz([peerdirec,presentfiles{anatfi}]);
        end
        weights=repmat(1.5,length(presentfiles),1);

        % add FA if present ? add to beginning since will use last entry to
        % converge
        if exist([subdirec,options.prefs.fa2anat],'file') && exist([peerdirec,options.prefs.fa2anat],'file')
            spfroms=[{ea_niigz([subdirec,options.prefs.fa2anat])},spfroms];
            sptos=[{ea_niigz([peerdirec,options.prefs.fa2anat])},sptos];

            weights=[0.5;weights];

        end

        defoutput=[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname];
        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep],'file')
            mkdir([subdirec,'MAGeT',filesep,'warps',filesep]);
        end
        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5'],'file')
            try
                ea_ants_nonlinear(sptos,spfroms,[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'.nii'],weights,options);
            catch
                ea_error(['Something went wrong - could not generate a nonlinear warp from ',subdirec,' to ',peerdirec,'.']);
            end
        end
        delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'.nii']); % we only need the warp

        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5'],'file')
            ea_error(['Something went wrong - could not generate a nonlinear warp from ',subdirec,' to ',peerdirec,'.']);
        end

        % Now export composite transform from MNI -> Peer -> Subject

        % check if peertransforms are available, if not, build them:
        if ~exist([peerfolders{peer},filesep,'glanatComposite',ea_getantstransformext([peerfolders{peer},filesep])],'file') || ...
                ~exist([peerfolders{peer},filesep,'glanatInverseComposite',ea_getantstransformext([peerfolders{peer},filesep])],'file')
            ea_normalize_ants(poptions,0);
        end

        antsdir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
        if ispc
            applyTransforms = ea_path_helper([antsdir, 'antsApplyTransforms.exe']);
        else
            applyTransforms = [antsdir, 'antsApplyTransforms.', computer('arch')];
        end

        template=ea_niigz([ea_space(options),'t2.nii']);
        prenii=ea_niigz([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
        cmd=[applyTransforms,' -r ',template,...
            ' -v ',...
            ' -t ',ea_path_helper([peerfolders{peer},filesep,'glanatComposite',ea_getantstransformext([peerfolders{peer},filesep])]),...
            ' -t ',ea_path_helper([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5']),...
            ' -o [',ea_path_helper([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2mni.nii']),',1]']; % temporary write out uncompressed (.nii) since will need to average slice by slice lateron.
        icmd=[applyTransforms,' -r ',prenii,...
            ' -v ',...
            ' -t ',ea_path_helper([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'InverseComposite.h5']),...
            ' -t ',ea_path_helper([peerfolders{peer},filesep,'glanatInverseComposite',ea_getantstransformext([peerfolders{peer},filesep])]),...
            ' -o [',ea_path_helper([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2sub.nii']),',1]']; % temporary write out uncompressed (.nii) since will need to average slice by slice lateron.
        if ~ispc
            system(['bash -c "', cmd, '"']);
            system(['bash -c "', icmd, '"']);
        else
            system(cmd);
            system(icmd);
        end

        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2mni.nii'],'file') || ...
                ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2sub.nii'],'file')
            ea_error(['Something went wrong - could not generate a nonlinear warp from ',subdirec,' to ',peerdirec,'.']);
        end
        % delete intermediary transforms
        delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5']);
        delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'InverseComposite.h5']);

    end

end
dosecondpass=1;

% aggregate all warps together and average them
clear ficell
warpbase=[options.root,options.patientname,filesep,'MAGeT',filesep,'warps',filesep];
delete([warpbase,'ave2mni.nii']);
fis=dir([warpbase,'*2mni.nii']);
for fi=1:length(fis)
    ficell{fi}=ea_niigz([warpbase,fis(fi).name]);
end
ea_robustaverage_nii(ficell,[warpbase,'ave2mni.nii']);
% need to do the following due to file format issues:
niione=ea_load_untouch_nii(ficell{fi});
niiave1=ea_load_untouch_nii([warpbase,'ave2mni_1.nii']);
niiave2=ea_load_untouch_nii([warpbase,'ave2mni_2.nii']);
niiave3=ea_load_untouch_nii([warpbase,'ave2mni_3.nii']);
niione.img(:,:,:,1,1)=niiave1.img;
niione.img(:,:,:,1,2)=niiave2.img;
niione.img(:,:,:,1,3)=niiave3.img;
ea_save_untouch_nii(niione,[warpbase,'ave2mni.nii']);

todelete=[];
if dosecondpass % discard warps that are too far off the robustmean
    for fi=1:length(ficell)
        testnii=ea_load_untouch_nii(ficell{fi});
        R=corr(testnii.img(1:100:end)',niione.img(1:100:end)');
        if R<0.9
            todelete=[todelete,fi];
        end
    end
    % redo the whole thing
    ficell(todelete)=[];
    ea_robustaverage_nii(ficell,[warpbase,'ave2mni.nii']);
    % need to do the following due to file format issues:
    niione=ea_load_untouch_nii(ficell{1});
    niiave1=ea_load_untouch_nii([warpbase,'ave2mni_1.nii']);
    niiave2=ea_load_untouch_nii([warpbase,'ave2mni_2.nii']);
    niiave3=ea_load_untouch_nii([warpbase,'ave2mni_3.nii']);
    niione.img(:,:,:,1,1)=niiave1.img;
    niione.img(:,:,:,1,2)=niiave2.img;
    niione.img(:,:,:,1,3)=niiave3.img;
    ea_save_untouch_nii(niione,[warpbase,'ave2mni.nii']);
end

gzip([warpbase,'ave2mni.nii']);

clear ficell
fis=dir([warpbase,'*sub.nii']);
for fi=1:length(fis)
    ficell{fi}=[warpbase,fis(fi).name];
end
ficell(todelete)=[];
ea_robustaverage_nii(ficell,[warpbase,'ave2sub.nii']);
% need to do the following due to file format issues:
niione=ea_load_untouch_nii(ficell{1});
niiave1=ea_load_untouch_nii([warpbase,'ave2sub_1.nii']);
niiave2=ea_load_untouch_nii([warpbase,'ave2sub_2.nii']);
niiave3=ea_load_untouch_nii([warpbase,'ave2sub_3.nii']);
niione.img(:,:,:,1,1)=niiave1.img;
niione.img(:,:,:,1,2)=niiave2.img;
niione.img(:,:,:,1,3)=niiave3.img;
ea_save_untouch_nii(niione,[warpbase,'ave2sub.nii']);
%
gzip([warpbase,'ave2sub.nii']);

movefile([warpbase,'ave2mni.nii.gz'],[subdirec,'glanatComposite.nii.gz']);
movefile([warpbase,'ave2sub.nii.gz'],[subdirec,'glanatInverseComposite.nii.gz']);

% delete older .h5 transforms if present.
ea_delete([subdirec,'glanatComposite.h5']);

ea_delete([subdirec,'glanatInverseComposite.h5']);

% apply warps as always:
ea_apply_normalization(options);


% finally, cleanup.

rmdir([subdirec,'MAGeT'],'s');


