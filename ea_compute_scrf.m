function ea_compute_scrf(handles)

options = getappdata(handles.scrf,'options');

% Hard-code the coregistration method to ANTs
options.coregmr.method = 'ANTs (Avants 2008)';

ea_prepareImages(options); % Redo each time since coreg may have changed.

if ~handles.mask0.Value % Use mask
    % Masks to be used
    masks = {options.subj.brainshift.anat.secondstepmask
        options.subj.brainshift.anat.thirdstepmask};
    if ~all(isfile(masks)) || options.overwriteapproved
        ea_createscrfmask(options);
    end
end

if handles.mask1.Value % Coarse mask
    masks = masks(1); % only use first mask.
end

if ~exist('masks','var')
    masks = {};
end

% Do subcortical coreigstration
affineTransform = ea_coregimages(options, options.subj.brainshift.anat.moving, ...
    options.subj.brainshift.anat.anchor, ...
    options.subj.brainshift.anat.scrf, ...
    {}, 1, masks);

% Rename transformation file
ea_mkdir(fileparts(options.subj.brainshift.transform.instore));
movefile(affineTransform{1}, options.subj.brainshift.transform.instore);
delete(affineTransform{2});

% Generated brainshift corrected postop images (coreg and norm)
ea_genscrfimages(options.subj);

% Refresh scrf status
ea_refreshscrf(options, handles);

% Dump method
ea_dumpmethod(options, 'brainshift');


% Prepare anchor and moving images for brain shift correction
function ea_prepareImages(options)
if ~isfile(options.subj.brainshift.anat.anchor) || options.overwriteapproved
    ea_mkdir([options.subj.brainshiftDir, filesep, 'anat']);

    from{1} = [ea_space,'bb.nii'];
    to{1} = [options.subj.brainshiftDir, filesep, 'anat', filesep, 'bb.nii'];

    % Warp template bb file into anchor native space
    try
        ea_apply_normalization_tofile(options,from,to,1,1);
    catch
        ea_error('Please perform normalization first.');
    end

    % Crop and reslice bb file
    ea_crop_nii(to{1});
    ea_reslice_nii(to{1}, to{1}, [0.4,0.4,0.4]);

    % Reslice anchor image for brainshift correction
    copyfile(options.subj.preopAnat.(options.subj.AnchorModality).coreg, options.subj.brainshift.anat.anchor)
    ea_conformspaceto(to{1}, options.subj.brainshift.anat.anchor);

    % Cleanup bb file
    delete(to{1});
end

% In case using CT, apply tone-mapping if needed, keep only one image
if strcmp(options.subj.postopModality, 'CT')
    if strcmp(options.prefs.scrf.tonemap, 'tp_')
        if (~isfile(options.subj.coreg.anat.postop.tonemapCT) || options.overwriteapproved) && isfile(options.subj.coreg.anat.postop.CT)
            ea_tonemapct(options, 'native');
        end
        options.subj.coreg.anat.postop = rmfield(options.subj.coreg.anat.postop, 'CT');
    else
        options.subj.coreg.anat.postop = rmfield(options.subj.coreg.anat.postop, 'tonemapCT');
    end
end

% Reslice post-op images for brainshift correction
postopCoregImage = struct2cell(options.subj.coreg.anat.postop);
postopMovingImage = strrep(postopCoregImage, options.subj.coregDir, options.subj.brainshiftDir);
for i = 1:length(postopCoregImage)
    if isfile(postopCoregImage{i})
        if ~isfile(postopMovingImage{i}) || options.overwriteapproved
            copyfile(postopCoregImage{i}, postopMovingImage{i});
            ea_conformspaceto(options.subj.brainshift.anat.anchor, postopMovingImage{i},1);
        end
    end
end

% In case post-op MRI, create mean image based on available images
if strcmp(options.subj.postopModality, 'MRI')
    if ~isfile(options.subj.brainshift.anat.moving) || options.overwriteapproved
        nii = cellfun(@(x) ea_load_nii(x), postopMovingImage);

        % Threshold
        for i=1:length(nii)
            nii(i).img(abs(nii(i).img)<0.1) = nan;
        end
    
        % Get mean image
        nii(1).img = ea_nanmean(cat(4, nii(:).img), 4);
        nii(1).img(isnan(nii(1).img)) = 0;
    
        % Save moving image
        nii(1).fname = options.subj.brainshift.anat.moving;
        ea_write_nii(nii(1));
    
        % Cleanup
        ea_delete(postopMovingImage);
    end
end
