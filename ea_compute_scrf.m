function ea_compute_scrf(handles)

options = getappdata(handles.scrf,'options');

% Hard-code the coregistration method to ANTs
options.coregmr.method = 'ANTs (Avants 2008)';

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

% Refresh scrf status
ea_refreshscrf(options, handles);
