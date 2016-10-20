function varargout=ea_normalize_maget(options)

if ischar(options) % return name of method.
    varargout{1}='MAGeT Brain-like Normalization (Chakravarty 2013)';
    varargout{2}={'SPM8','SPM12'};
    varargout{3}=1; % hassettings.
    return
end

reforce=1;

peerfolders=ea_getmagetpeers(options);


%% step 0: check if all subjects have been processed with an ANTs-based normalization function
for peer=1:length(peerfolders)
    poptions=options;
    [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
    poptions.root=[poptions.root,filesep];
    poptions=ea_assignpretra(poptions);
    % make sure peer has been normalized using ANTs
    if ~ismember(ea_whichnormmethod([peerfolders{peer},filesep]),ea_getantsnormfuns)
        ea_dumpnormmethod(poptions,'ea_normalize_ants_multimodal');
        ea_normalize_ants_multimodal(poptions)
    else
        % make sure peer's anatomy files have been coregistered.
        if ~ea_seemscoregistered(poptions)
            ea_coreg_all_mri(poptions,0);
        end
    end
end

subdirec=[options.root,options.patientname,filesep];
if ~ea_seemscoregistered(options)
    ea_coreg_all_mri(options,0);
end

%% step 1, setup DISTAL warps back to sub via each peer brain
earoot=ea_getearoot;

for peer=1:length(peerfolders)
    
    clear spfroms sptos weights metrics
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
        if ~isequal(subpresentfiles,peerpresentfiles) % then I made something wrong.
            keyboard
        else
            presentfiles=subpresentfiles;
            clear subpresentfiles peerpresentfiles
        end
        clear sptos spfroms metrics weights
        for anatfi=1:length(presentfiles)
            spfroms{anatfi}=ea_niigz([subdirec,presentfiles{anatfi}]);
            sptos{anatfi}=ea_niigz([peerdirec,presentfiles{anatfi}]);
            metrics{anatfi}='MI';
        end
        weights=repmat(1.5,length(presentfiles),1);
        
        % add FA if present ? add to beginning since will use last entry to
        % converge
        if exist([subdirec,options.prefs.fa2anat],'file') && exist([peerdirec,options.prefs.fa2anat],'file')
            spfroms=[{ea_niigz([subdirec,options.prefs.fa2anat])},spfroms];
            sptos=[{ea_niigz([peerdirec,options.prefs.fa2anat])},sptos];
            
            weights=[0.5;weights];
            metrics=[{'MI'},metrics];
            
        end
        
        defoutput=[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname];
        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep],'file')
            mkdir([subdirec,'MAGeT',filesep,'warps',filesep]);
        end
        
        try
            ea_ants_nonlinear(sptos,spfroms,[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'.nii'],weights,metrics,options);
        catch
            ea_error(['Something went wrong - could not generate a nonlinear warp from ',subdirec,' to ',peerdirec,'.']);
        end
        delete([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'.nii']); % we only need the warp
        
        if ~exist([subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5'],'file')
            ea_error(['Something went wrong - could not generate a nonlinear warp from ',subdirec,' to ',peerdirec,'.']);
        end
        
        % Now export composite transform from MNI -> Peer -> Subject
        
        if ispc
            sufx='.exe';
        else
            sufx=computer('arch');
        end
        
        antsApply=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep,'antsApplyTransforms.',sufx];
        
        template=ea_niigz([ea_getearoot,'templates',filesep,'mni_hires.nii']);
        prenii=ea_niigz([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
        cmd=[antsApply,' -r ',template,' -t ',[peerfolders{peer},filesep,'glanatComposite',ea_getantstransformext([peerfolders{peer},filesep],options)],' -t ',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'Composite.h5'],' -o [',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2mni.nii',',1]']]; % temporary write out uncompressed (.nii) since will need to average slice by slice lateron.
        icmd=[antsApply,' -r ',prenii,' -t ',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'InverseComposite.h5'],' -t ',[peerfolders{peer},filesep,'glanatInverseComposite',ea_getantstransformext([peerfolders{peer},filesep],options)],' -o [',[subdirec,'MAGeT',filesep,'warps',filesep,poptions.patientname,'2sub.nii',',1]']]; % temporary write out uncompressed (.nii) since will need to average slice by slice lateron.
        if ~ispc
            system(['bash -c "', cmd, '"']);
            system(['bash -c "', icmd, '"']);
        else
            system(cmd);
            system(icmd);
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
if exist([subdirec,'glanatComposite.h5'],'file')
    delete([subdirec,'glanatComposite.h5']);
end

if exist([subdirec,'glanatInverseComposite.h5'],'file')
    delete([subdirec,'glanatInverseComposite.h5']);
end

% % now convert to .h5 again and place in sub directory:
% antsApply=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep,'antsApplyTransforms.',sufx];
% template=[ea_getearoot,'templates',filesep,'mni_hires.nii'];
% prenii=[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized];
% cmd=[antsApply,' -r ',template,' -t ',[warpbase,'ave2mni.nii.gz'],' -o [',[subdirec,'glanatComposite.nii.gz,1]']];
% icmd=[antsApply,' -r ',prenii,' -t ',[warpbase,'ave2sub.nii.gz'],' -o [',[subdirec,'glanatInverseComposite.nii.gz,1]']];
% if ~ispc
%     system(['bash -c "', cmd, '"']);
%     system(['bash -c "', icmd, '"']);
% else
%     system(cmd);
%     system(icmd);
% end

% apply warps as always:
ea_apply_normalization(options);


% finally, cleanup.

rmdir([subdirec,'MAGeT'],'s');



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

