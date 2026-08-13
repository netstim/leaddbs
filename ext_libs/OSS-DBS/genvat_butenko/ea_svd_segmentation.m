function ea_svd_segmentation(images, anchorImage, segmenter, segmaskPath, alg4PVS,alg4WMH,alg4Lacunes,sub_ID,only_segmask,alg4MB)
% Add enlarged perivascular spaces, WMH and lacunes to segmask
% By Butenko, konstantinmgtu@gmail.com
arguments
    images                          % struct of paths to images
    anchorImage                     % path to the image which is used as anchor
    %coregister = true               % if true, all images co-registered to anchorImage
    segmenter = 'SynthSeg'          % algorithm for general segmentation
    segmaskPath = 'default'         % segmask output path    %
    alg4PVS = 'segcsvd'             % algorithm for PVS segmentation
    alg4WMH = 'segcsvd'             % algorithm for WMH segmentation
    alg4Lacunes = 'TeamTea'         % algorithm for lacune segmentation
    sub_ID = 'JD'                   % subject label
    only_segmask = false            % only keep segmask.nii, remove all intermediate files
    alg4MB = 'None'
end

% for now, enforce co-registration
coregister = false;
%anchorNii = anchorImage;

% check if all images provided actually exist
for im_i = 1:size(images,1)
    image = [images(im_i).folder,filesep,images(im_i).name];
    if ~isfile(image)
        fprintf('\nFor %s, not all listed modalities were found\n',sub_ID)
        return
    end
end

% logic check
if strcmp(alg4PVS,'segcsvd') & ~strcmp(alg4WMH,'segcsvd')
    disp("When using segcsvd for PVS, also use it for WMH")
    return
end

% autodefinitions
if strcmp(segmaskPath,'default')
    [workingDir,anchorName,~] = fileparts(anchorImage); 
    if strcmp(anchorImage(end-2:end),'.gz')
        anchorName = anchorName(1:end-4);
    end
    segmaskPath = [workingDir, filesep, 'segmask.nii'];
else
    [workingDir,~,~] = fileparts(segmaskPath);
    [~,anchorName,~] = fileparts(anchorImage); 
end

% clean-up
ea_delete([workingDir,filesep,'anchor_synthSeg.nii.gz']);
ea_delete([workingDir,filesep,'T1_synthSeg.nii.gz']);
ea_delete([workingDir,filesep,'segmask_synthSeg.nii']);

segmaskFile = [workingDir,filesep,'segmask_',segmenter,'.nii'];
segmaskSVDFile = segmaskPath;

% we need .nii for spm commands
if strcmp(anchorImage(end-2:end),'.gz')
    % SPM does not handle .gz
    gunzip(anchorImage);
    anchorNii = [anchorImage(1:end-3)];
    % for clean-up
    images2remove{1,1} = anchorNii;
    cnt_rem = 2;
else
    anchorNii = anchorImage;
end

% 0) additional co-registration
% also re-slices to the same voxel space, which is handy
if coregister
    for im_i = 1:size(images,1)
        image = [images(im_i).folder,filesep,images(im_i).name];

        if ~strcmp(image,anchorImage)
            if strcmp(image(end-2:end),'.gz')
                % SPM does not handle .gz
                gunzip(image);
                image = [image(1:end-3)];
                images2remove{cnt_rem,1} = image;
                cnt_rem = cnt_rem + 1;
            end

            if ~strcmp(image,anchorImage)
                [~,imageName,imageFormat] = fileparts(image); 
    
                where_to_save_cormat = [workingDir,filesep,'spm_affine_',imageName,imageFormat];
            end
            ea_spm_coreg(where_to_save_cormat,image,anchorNii,'nmi',1, '',1);
            images(im_i).name = ['r',imageName,'.nii'];
        else
            images(im_i).name =  images(im_i).name(1:end-3);
        end
    end
end

if strcmp(segmenter,'SynthSeg')
    for im_i = 1:size(images,1)
        image2segment = [images(im_i).folder,filesep,images(im_i).name];
        [mask_nifti,synth_fn] = get_SynthSeg_segmask(workingDir,image2segment,anchorNii,segmaskFile);
        if nnz(mask_nifti.img) ~= 0
            % segmentation was "successful"
            break
        elseif im_i == size(images,1)
            disp("SynthSeg failed for all modalities, skipping the subject")
            return
        end
    end    
else
    disp("Only SynthSeg is currently supported for segmentation")
    return
end

% check which modalities we have, mask to brain
found = [false;false;false];  % we might have T1, T2 and FLAIR
%mask_nifti = ea_load_nii(segmaskFile);

% Store all co-registed and brain masked files in one directory
TT_input_path = [workingDir,filesep,'TT_input'];
TT_output_path = [workingDir,filesep,'TT_output'];
ea_delete(TT_input_path);
%ea_delete(TT_output_path);
mkdir(TT_input_path);
mkdir(TT_output_path);

for im_i = 1:size(images,1)
    % we need T1, T2 and FLAIR
    image = [images(im_i).folder,filesep,images(im_i).name];

    % check if subject prefix is already used in the naming
    %subj_prefix = ['sub-',sub_ID];
    % if contains(sub_ID,'rsub-')
    %     subj_prefix = ['sub-',sub_ID(6:end)];
    if contains(sub_ID,'sub-')
        subj_prefix = sub_ID;
    else
        subj_prefix = ['sub-',sub_ID];
    end

    % apply the mask
    image_nifti = ea_load_nii(image);
    image_nifti.img(mask_nifti.img < 0.5) = 0.0;

    if coregister
        ea_delete(image)  % comment this out if we need .nii somewhere downstream  
    end
    images(im_i).folder = TT_input_path;

    % rename to TeamTea notation (synthSeg, segmsvd and MARS support it directly)
    if contains(image, 'T1w') || contains(image, 'anat_t1') || contains(image, 'masked_T1') 
        image_T1 = [subj_prefix,'_space-anchor_desc-masked_T1.nii.gz'];
        images(im_i).name = image_T1;
        found(1) = true;
    elseif contains(image, 'FLAIR') || contains(image, 'anat_flair') 
        image_FLAIR = [subj_prefix,'_space-anchor_desc-masked_FLAIR.nii.gz'];
        images(im_i).name = image_FLAIR;
        found(3) = true;
    elseif contains(image, 'T2starw') || contains(image, 'T2wS') || contains(image, 'anat_t2S') || contains(image, 'masked_T2S')
        image_T2S = [subj_prefix,'_space-anchor_desc-masked_T2S.nii.gz'];
        images(im_i).name = image_T2S;
    elseif contains(image, 'T2w') || contains(image, 'anat_t2') || contains(image, 'masked_T2')
        image_T2 = [subj_prefix,'_space-anchor_desc-masked_T2.nii.gz'];
        images(im_i).name = image_T2;
        found(2) = true;
    else
        % we do not need this image
        continue
    end

    image_nifti.fname = [images(im_i).folder,filesep,images(im_i).name];
    ea_write_nii(image_nifti);

end


if strcmp(alg4Lacunes,'TeamTea')

    if ~all(found)
        fprintf('\nFor %s, not all required modalities were found for TeamTea Lacune segmentation\n',sub_ID)
        return
    end

    % D = gpuDevice;
    % reset(D);

    disp("Running lacune segmentation with TeamTea")

    % % run TeamTea docker
    % you might need to add user to the docker group

    system(['docker run --gpus all --rm ', ...
        '-v ', TT_input_path, ':/input ', ...
        '-v ', TT_output_path, ':/output ', ...
        '59f961f4a16d']);    % TeamTea
    
    % for now run externally
    %fprintf('\nIn the terminal: sudo docker run --rm --gpus all -v %s:/input -v %s:/output 59f961f4a16d\n',TT_input_path,TT_output_path)
    %disp("Copy the docker command and make sure the process is over")

    lacunesFile = [TT_output_path,filesep,'images',filesep,subj_prefix,'_space-anchor_desc-prediction.nii.gz'];
    if ~isfile(lacunesFile)
        disp("Lacune segmentation with TeamTea failed!")
        disp(sub_ID)
    else
        lacunesFileAccessible = [workingDir,filesep,subj_prefix,'_space-anchor_desc-lacunes.nii.gz'];
        copyfile(lacunesFile,lacunesFileAccessible)
    end

    % clean-up
    %ea_delete(TT_input_path);
else
    disp("alg4Lacunes options is not recognized")
    %return
end

if strcmp(alg4PVS,'TeamTea')
    if ~all(found)
        fprintf('\nFor %s, not all required modalities were found for TeamTea PVS segmentation\n',sub_ID)
        return
    end

    disp("Running PVS segmentation with TeamTea")
    % % run TeamTea docker
    % you might need to add user to the docker group

    system(['docker run --gpus all --rm ', ...
        '-v ', TT_input_path, ':/input ', ...
        '-v ', TT_output_path, ':/output ', ...
        '13387316fdfe']);    % TeamTea    

    pvsFile = [TT_output_path,filesep,'images',filesep,subj_prefix,'_space-anchor_desc-prediction.nii.gz'];
    if ~isfile(pvsFile)
        disp("PVS segmentation with TeamTea failed!")
        disp(sub_ID)
    else
        pvsFileAccessible = [workingDir,filesep,subj_prefix,'_space-anchor_desc-PVS.nii.gz'];
        copyfile(pvsFile,pvsFileAccessible)
    end
end

if strcmp(alg4MB,'TeamTea')
    if ~all(found)
        fprintf('\nFor %s, not all required modalities were found for TeamTea microbleed segmentation\n',sub_ID)
        return
    end

    % % pretend that T2* is T2
    % image_T2S = [subj_prefix,'_space-anchor_desc-masked_T2S.nii.gz'];
    % copyfile([TT_input_path,filesep,image_T2],[TT_input_path,filesep,image_T2S])

    % % run TeamTea docker
    % you might need to add user to the docker group
    disp("Running microbleed segmentation with TeamTea")

    system(['docker run --gpus all --rm ', ...
        '-v ', TT_input_path, ':/input ', ...
        '-v ', TT_output_path, ':/output ', ...
        '95a63079481e']);    % TeamTea    

    mbFile = [TT_output_path,filesep,'images',filesep,subj_prefix,'_space-anchor_desc-prediction.nii.gz'];
    if ~isfile(mbFile)
        disp("MB segmentation with TeamTea failed!")
        disp(sub_ID)
    else
        mbFileAccessible = [workingDir,filesep,subj_prefix,'_space-anchor_desc-MB.nii.gz'];
        copyfile(mbFile,mbFileAccessible)
    end
end

if strcmp(segmenter,'SynthSeg')
    % keep all in one place
    copyfile([workingDir,filesep,synth_fn],TT_input_path)
end

% 2) do segmsvd WMH and PVS
if strcmp(alg4WMH,'segcsvd')

    if ~found(3)
        disp("No FLAIR image for WMH segmentation with segcsvd !")
        return
    end

    disp("Running WMH segmentation with segcsvd")

    segcsvdDir = [workingDir,filesep,'segcsvd'];
    mkdir(segcsvdDir);

    inputDir = [workingDir,filesep,'TT_input'];

    % synthSeg definitions
	%synth_fn = [anchorName,'_synthSeg.nii.gz'];
	seg_wmh_fn = 'seg_wmh.nii.gz';
    seg_wmh_thr_float = 0.1;
	seg_wmh_thr = char(string(seg_wmh_thr_float));
	skip_mask_and_bias = 'false';
	cleanup = 'true';
	
	%rm -f /${out_dir}/${seg_wmh_fn}
    ea_delete([segcsvdDir,filesep,seg_wmh_fn])
     
    system(['docker run --gpus all --rm ', ...
        '-v ',  inputDir, ':/indir ', ...
        '-v ', segcsvdDir, ':/outdir ', ...
        '-w / segcsvd_rc03 segment_wmh ', ...
        '/indir/',image_FLAIR, ...
        ' /indir/',synth_fn, ...
        ' /outdir/',seg_wmh_fn, ...
        ' 1 "96,128" ', seg_wmh_thr, ' 1 ', ...
        skip_mask_and_bias, ' ', cleanup]);

    %fprintf('\nIn the terminal: sudo docker run --rm --gpus all -v %s:/indir -v %s:/outdir -w / segcsvd_rc03 segment_wmh /indir/%s /indir/%s /outdir/%s 1 "96,128" %s 1 %s %s\n',inputDir,segcsvdDir,image_FLAIR,synth_fn,seg_wmh_fn,seg_wmh_thr,skip_mask_and_bias,cleanup)
	%disp("Copy the docker command and make sure the process is over")

    if ~isfile([segcsvdDir,filesep,seg_wmh_fn])
        disp("WMH segmentation with segcsvd failed!")
        disp(sub_ID)
    else
        wmhFile = [segcsvdDir,filesep,seg_wmh_fn];
    end

    if strcmp(alg4PVS,'segcsvd')

        if ~found(1)
            disp("No T1 image for PVS segmentation with segcsvd !")
            %return
        end

        disp("Running PVS segmentation with segcsvd")

        seg_pvs_thr_float = 0.35;
	    seg_pvs_thr = char(string(seg_pvs_thr_float));
	    seg_pvs_fn = 'seg_pvs.nii.gz';

        ea_delete([segcsvdDir,filesep,seg_pvs_fn])

        system(['docker run --gpus all --rm ', ...
            '-v ',  inputDir, ':/indir ', ...
            '-v ', segcsvdDir, ':/outdir ', ...
            '-w / segcsvd_rc03 segment_pvs ', ...
            '/indir/',image_T1, ...
            ' /indir/',synth_fn, ...
            ' /outdir/',seg_wmh_fn, ...
            ' /outdir/',seg_pvs_fn, ...
            ' "1.0,1.4" 0 ', 'true', ' ', ...
            cleanup, ' ',seg_pvs_thr]);

        %fprintf('\nIn the terminal: sudo docker run --rm --gpus all -v %s:/indir -v %s:/outdir -w / segcsvd_rc03 segment_pvs /indir/%s /indir/%s /outdir/%s /outdir/%s "1.0,1.4" 0 %s %s %s\n',inputDir,segcsvdDir,image_T1,synth_fn,seg_wmh_fn,seg_pvs_fn,skip_mask_and_bias,cleanup,seg_pvs_thr)
    	%disp("Copy the docker command and make sure the process is over")

        pvsFileAccessible = [segcsvdDir,filesep,seg_pvs_fn];
        if ~isfile(pvsFileAccessible)
            disp("PVS segmentation with segcsvd failed!")
            disp(sub_ID)
        end
    end

end


if strcmp(alg4PVS,'PINGU')
    % TBA
    % DataRaw and _0000.nii.gz
end


% burn in lacunes and PVS as CSF
segmask = ea_load_nii(segmaskFile);

% WMH, PVS and Lacune burn-in is ordered based on the VCM relevance
if ~strcmp(alg4WMH,'None')
    try
        wmh = ea_load_nii(wmhFile);
        if all(size(segmask.img) == size(wmh.img)) & all(abs(segmask.mat - wmh.mat) <= 0.001, 'all')
            % WMH ca be burned in as default to be filtered out during PAM
            segmask.img(wmh.img >= seg_wmh_thr_float) = 0;
        else
            disp("Voxel space of segmask and WMH did not match, check!")
        end
    catch ME
        fprintf('\nNo WMH segmentation found for %s\n',sub_ID)
    end
end

if ~strcmp(alg4PVS,'None')
    try
        pvs = ea_load_nii(pvsFileAccessible);
        if all(size(segmask.img) == size(pvs.img)) & all(abs(segmask.mat - pvs.mat) <= 0.001, 'all')
            if strcmp(alg4PVS,'segcsvd')
                segmask.img(pvs.img >= seg_pvs_thr_float) = 3;
            elseif strcmp(alg4PVS,'TeamTea')
                segmask.img(pvs.img == 1) = 3;
            end
        else
            disp("Voxel space of segmask and PVS did not match, check!")
        end
    catch ME
        fprintf('\nNo PVS segmentation found for %s\n',sub_ID)
    end
end

if ~strcmp(alg4Lacunes,'None')
    % file cannot be loaded if no lacunes found?
    try
        lacunes = ea_load_nii(lacunesFileAccessible);
        if all(size(segmask.img) == size(lacunes.img)) & all(abs(segmask.mat - lacunes.mat) <= 0.001, 'all')
            segmask.img(lacunes.img >= 0.99) = 3;
        else
            disp("Voxel space of segmask and lacunes did not match, check!")
        end
    catch ME
        fprintf('\nNo lacune segmentation found for %s \n',sub_ID)
    end
end

% burn in microbleeds as default (not conductive since hypointense on T2?)
if ~strcmp(alg4MB,'None')
    try
        microbleeds = ea_load_nii(mbFileAccessible);
        if all(size(segmask.img) == size(microbleeds.img)) & all(abs(segmask.mat - microbleeds.mat) <= 0.001, 'all')
            segmask.img(microbleeds.img >= 0.99) = 0;
        else
            disp("Voxel space of segmask and microbleeds did not match, check!")
        end
    catch ME
        fprintf('\nNo microbleed segmentation found for %s \n',sub_ID)
    end
end


segmask.fname = segmaskSVDFile;
ea_write_nii(segmask)

function [segmask_nifti,synth_fn] = get_SynthSeg_segmask(workingDir,image2segment,anchorNii,segmaskFile)

    %synth_fn = 'T1_synthSeg.nii.gz';
    %synth_fn = 'anchor_SynthSeg.nii.gz';

    % it is in .nii after co-registration
    [~,image2segment_name,~] = fileparts(image2segment);
    synth_fn = [image2segment_name,'_SynthSeg.nii.gz'];

    multisegment = [workingDir,filesep,synth_fn];  % output of the segmentation alg.

    % if strcmp(alg4Lacunes,'TeamTea')
    %     T12segment = [TT_input_path,filesep,image_T1];
    % else
    %     T12segment = [workingDir,filesep,image_T1];
    % end
    % ea_synthseg(T12segment, multisegment)

    % if anchor fails, try another image
    if ~isfile(multisegment)
        ea_synthseg(image2segment, multisegment) 
    end
    ea_convert_synthSeg2segmask(multisegment, segmaskFile);

    % reslice segmask to AnchorImage space
    % does not matter which image, they are all co-registered
    % but maybe better use anchor in the 

    % gunzip(T12segment)
    % T12segment_nii = T12segment(1:end-3);
    % ea_conformspaceto(T12segment_nii,segmaskFile)

    % if strcmp(anchorImage(end-2:end),'.gz') & ~isfile(anchorNii)
    %     % SPM does not handle .gz
    %     gunzip(anchorImage);
    %     anchorNii = [anchorImage(1:end-3)];
    % end
    ea_conformspaceto(anchorNii,segmaskFile,0)
    segmask_nifti = ea_load_nii(segmaskFile);
end

% always clean-up gunzipped files (if they were gzipped originally)
if strcmp(anchorImage(end-2:end),'.gz') && coregister
    for i = 1:size(images2remove,1)
        ea_delete(images2remove{i,1});
    end
end

if only_segmask
    ea_delete(segmaskFile)
    ea_delete([workingDir,filesep,synth_fn])
    ea_delete(TT_input_path)
    %ea_delete(TT_output_path)
    ea_delete(lacunesFileAccessible)
    ea_delete(pvsFileAccessible)
    ea_delete(wmhFile)
    ea_delete(segcsvdDir)
    ea_delete([workingDir,filesep,'spm_affine_*'])

    % if strcmp(anchorImage(end-2:end),'.gz')
    %     ea_delete(anchorNii);
    % end
end

end
