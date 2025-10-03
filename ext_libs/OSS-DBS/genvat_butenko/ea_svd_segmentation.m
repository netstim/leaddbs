function ea_svd_segmentation(images, anchorImage, coregister, segmenter, segmaskPath, alg4PVS,alg4WMH,alg4Lacunes,sub_ID)
% Add enlarged perivascular spaces, WMH and lacunes to segmask
% By Butenko, konstantinmgtu@gmail.com
arguments
    images                          % struct of paths to images
    anchorImage                     % path to the image which is used as anchor
    coregister = true               % if true, all images co-registered to anchorImage
    segmenter = 'SynthSeg'          % algorithm for general segmentation
    segmaskPath = 'default'         % segmask output path    %
    alg4PVS = 'segcsvd'             % algorithm for PVS segmentation
    alg4WMH = 'segcsvd'             % algorithm for WMH segmentation
    alg4Lacunes = 'TeamTea'         % algorithm for lacune segmentation
    sub_ID = 'JD'                   % subject label
end

% better use T2 as anchor to avoid segcsvd PVS bug

% logic check
if strcmp(alg4PVS,'segcsvd') & ~strcmp(alg4WMH,'segcsvd')
    ea_warndlg("When using segcsvd for PVS, also use it for WMH")
    return
end

% autodefinitions
if strcmp(segmaskPath,'default')
    
    [workingDir,anchorName,~] = fileparts(anchorImage); 
    if strcmp(anchorImage(end-2:end),'.gz')
        anchorName = anchorName(1:end-4);
    end
else
    [workingDir,~,~] = fileparts(segmaskPath);
    [~,anchorName,~] = fileparts(anchorImage); 
end

segmaskFile = [workingDir,filesep,'segmask_synthSeg.nii'];
segmaskSVDFile = [workingDir,filesep,'segmask.nii'];

% 0) additional co-registration
% also re-slices to the same voxel space, which is handy
if coregister

    if strcmp(anchorImage(end-2:end),'.gz')
        % SPM does not handle .gz
        gunzip(anchorImage);
        anchorNii = [anchorImage(1:end-3)];
    else
        anchorNii = anchorImage;
    end

    for im_i = 1:size(images,1)
        image = [images(im_i).folder,filesep,images(im_i).name];

        if ~strcmp(image,anchorImage)
            if strcmp(image(end-2:end),'.gz')
                % SPM does not handle .gz
                gunzip(image);
                image = [image(1:end-3)];
            end

            if ~strcmp(image,anchorImage)
                [~,imageName,imageFormat] = fileparts(image); 
    
                where_to_save_cormat = [workingDir,filesep,'spm_affine_',imageName,imageFormat];
            end
            ea_spm_coreg(where_to_save_cormat,image,anchorNii,'nmi',1, '',1)
            images(im_i).name = ['r',imageName,'.nii'];
        end
    end
end

found = [false;false;false];  % we might T1, T2 and FLAIR
if strcmp(alg4Lacunes,'TeamTea')
    % rename to TeamTea notation (synthSeg, segmsvd and MARS support it directly)
    TT_input_path = [workingDir,filesep,'TT_input'];
    TT_output_path = [workingDir,filesep,'TT_output'];

    mkdir(TT_input_path);
    mkdir(TT_output_path);

    for im_i = 1:size(images,1)
        % we need T1, T2 and FLAIR
        image = [images(im_i).folder,filesep,images(im_i).name];

        % images should be gzipped
        if ~strcmp(image(end-2:end),'.gz')
            gzip(image)
            ea_delete(image)  % comment this out if we need .nii somewhere downstream
            image = [image,'.gz'];

            % use gzip further
            images(im_i).name = [images(im_i).name,'.gz'];
        end
            
        if contains(image, 'T1w') || contains(image, 'anat_t1')
            image_T1 = ['sub-',sub_ID,'_space-default_desc-masked_T1.nii.gz'];
            copyfile(image,[TT_input_path,filesep,image_T1])
            found(1) = true;
        elseif contains(image, 'T2w') || contains(image, 'anat_t2')
            image_T2 = ['sub-',sub_ID,'_space-default_desc-masked_T2.nii.gz'];
            copyfile(image,[TT_input_path,filesep,image_T2])
            found(2) = true;
        elseif contains(image, 'FLAIR') || contains(image, 'anat_flair')
            image_FLAIR = ['sub-',sub_ID,'_space-default_desc-masked_FLAIR.nii.gz'];
            copyfile(image,[TT_input_path,filesep,image_FLAIR])
            found(3) = true;
        else
            % we do not need this image
            continue
        end
    end

    if ~all(found)
        fprintf('\nFor %s, not all required modalities were found for TeamTea Lacune segmentation\n',sub_ID)
        %return
    end

    % % run TeamTea docker
    % you might need to add user to the docker group
    system(['docker run --gpus all --rm ', ...
        '-v ', TT_input_path, ':/input ', ...
        '-v ', TT_output_path, ':/output ', ...
        '59f961f4a16d']);
    
    % for now run externally
    %fprintf('\nIn the terminal: sudo docker run --rm --gpus all -v %s:/input -v %s:/output 59f961f4a16d\n',TT_input_path,TT_output_path)
    %ea_warndlg("Copy the docker command and make sure the process is over")

    lacunesFile = [TT_output_path,filesep,'images',filesep,'sub-',sub_ID,'_space-default_desc-prediction.nii.gz'];
    if ~isfile(lacunesFile)
        ea_warndlg("Lacune segmentation with segcsvd failed!")
        disp(sub_ID)
    end

    % clean-up
    %ea_delete(TT_input_path);

elseif strcmp(alg4Lacunes,'None')
    % just classify images
    for im_i = 1:size(images,1)
        if contains(image, 'T1w') || contains(image, 'anat_t1')
            image_T1 = images(im_i).name;
            found(1) = true;
        elseif contains(image, 'T2w') || contains(image, 'anat_t2')
            image_T2 = images(im_i).name;
            found(2) = true;
        elseif contains(image, 'FLAIR') || contains(image, 'anat_flair')
            image_FLAIR = images(im_i).name;
            found(3) = true;
        else
            % we do not need this image
            continue
        end
    end
else
    ea_warndlg("alg4Lacunes options is not recognized")
    %return
end

% we can do it in a batch
% we can create a batch bash script from matlab
%create_bash_batch()

% 1) do synthseg (it is good with general anatomy)
% segment the anchor modality

if strcmp(segmenter,'SynthSeg')
    synth_fn = 'T1_synthSeg.nii.gz';
    multisegment = [workingDir,filesep,synth_fn];  % output of the segmentation alg.

    if ~isfile(multisegment)
        if strcmp(alg4Lacunes,'TeamTea')
            T12segment = [TT_input_path,filesep,image_T1];
        else
            T12segment = [workingDir,filesep,image_T1];
        end
        ea_synthseg(T12segment, multisegment)
        ea_convert_synthSeg2segmask(multisegment, segmaskFile);

        % reslice segmask to AnchorImage space
        % does not matter which image, they are all co-registered
        % but maybe better use anchor in the 
        gunzip(T12segment)
        T12segment_nii = T12segment(1:end-3);
        ea_conformspaceto(T12segment_nii,segmaskFile)

    end

    if strcmp(alg4Lacunes,'TeamTea')
        % keep all in one place
        copyfile(multisegment,[workingDir,filesep,'TT_input'])
    end
else
    ea_warndlg("Only SynthSeg is currently supported for segmentation")
    return
end

% 2) do segmsvd WMH and PVS
if strcmp(alg4WMH,'segcsvd')

    if ~found(3)
        ea_warndlg("No FLAIR image for WMH segmentation with segcsvd !")
        %return
    end

    segcsvdDir = [workingDir,filesep,'segcsvd'];
    mkdir(segcsvdDir);

    if strcmp(alg4Lacunes,'TeamTea')
        inputDir = [workingDir,filesep,'TT_input'];
    else
        inputDir = workingDir;
    end

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
	%ea_warndlg("Copy the docker command and make sure the process is over")

    if ~isfile([segcsvdDir,filesep,seg_wmh_fn])
        ea_warndlg("WMH segmentation with segcsvd failed!")
        disp(sub_ID)
    else
        wmhFile = [segcsvdDir,filesep,seg_wmh_fn];
    end

    if strcmp(alg4PVS,'segcsvd')

        if ~found(1)
            ea_warndlg("No T1 image for PVS segmentation with segcsvd !")
            %return
        end

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
            ' "1.0,1.4" 0 ', skip_mask_and_bias, ' ', ...
            cleanup, ' ',seg_pvs_thr]);

        %fprintf('\nIn the terminal: sudo docker run --rm --gpus all -v %s:/indir -v %s:/outdir -w / segcsvd_rc03 segment_pvs /indir/%s /indir/%s /outdir/%s /outdir/%s "1.0,1.4" 0 %s %s %s\n',inputDir,segcsvdDir,image_T1,synth_fn,seg_wmh_fn,seg_pvs_fn,skip_mask_and_bias,cleanup,seg_pvs_thr)
    	%ea_warndlg("Copy the docker command and make sure the process is over")

        pvsFile = [segcsvdDir,filesep,seg_pvs_fn];
        if ~isfile(pvsFile)
            ea_warndlg("PVS segmentation with segcsvd failed!")
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
        if all(size(segmask.img) == size(wmh.img)) & all(segmask.mat == wmh.mat)
            % WMH ca be burned in as default to be filtered out during PAM
            segmask.img(wmh.img >= seg_wmh_thr_float) = 0;
        else
            ea_warndlg("Voxel space of segmask and WMH did not match, check!")
        end
    catch ME
        fprintf('\nNo WMH segmentation found for %s\n',sub_ID)
    end
end

if ~strcmp(alg4PVS,'None')
    try
        pvs = ea_load_nii(pvsFile);
        if all(size(segmask.img) == size(pvs.img)) & all(segmask.mat == pvs.mat)
            segmask.img(pvs.img >= seg_pvs_thr_float) = 3;
        else
            ea_warndlg("Voxel space of segmask and PVS did not match, check!")
        end
    catch ME
        fprintf('\nNo PVS segmentation found for %s\n',sub_ID)
    end
end

if ~strcmp(alg4Lacunes,'None')
    % file cannot be loaded if no lacunes found?
    try
        lacunes = ea_load_nii(lacunesFile);
        if all(size(segmask.img) == size(lacunes.img)) && all(segmask.mat == lacunes.mat)
            segmask.img(lacunes.img >= 0.99) = 3;
        else
            ea_warndlg("Voxel space of segmask and lacunes did not match, check!")
        end
    catch ME
        fprintf('\nNo lacune segmentation found for %s \n',sub_ID)
    end
end

segmask.fname = segmaskSVDFile;
ea_write_nii(segmask)
