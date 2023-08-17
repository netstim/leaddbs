function ea_show_normalization(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ~isfield(options, 'leadprod')
    callingfunction = 'normalization dbs';
else
    callingfunction = ['normalization ', options.leadprod];
end

disp('Preparing images to show Normalization...');

try
    % Normalized anchor image
    checkImge = options.subj.norm.anat.preop.(options.subj.AnchorModality);

    % Set title string
    modality = ea_getmodality(checkImge);
    title = ['MNI ', upper(options.primarytemplate), ' (wireframes) & Preoperative MRI (',modality,')'];

    % Load and adjust wire frame
    load([ea_space, 'wires.mat'], 'wires');
    wires = single(wires);
    wires = wires/255;
    wires = wires .* 0.2;
    wires = wires + 0.8;

    % Load and adjust normalized anchor image
    subj = ea_load_nii(checkImge);
    subjImg = subj.img;
    subjImg = (subjImg-min(subjImg(:)))/(max(subjImg(:)));
    subjImg(subjImg>0.5) = 0.5;
    subjImg = (subjImg-min(subjImg(:)))/(max(subjImg(:)));

    % Load and adjust template image
    template = ea_load_nii([ea_space, options.primarytemplate, '.nii']);
    templateImg = template.img;
    templateImg(:) = zscore(templateImg(:));
    templateImg = (templateImg-min(templateImg(:)))/(max(templateImg(:))-min(templateImg(:)));

    % Set joint image
    jointImg = subjImg .* wires;
    jointImg = (jointImg - min(jointImg(:)))./(max(jointImg(:))-std(jointImg(:)));
    jointImg(jointImg>1) = 1;

    % Set data for detail viewer
    subjImg = single(subjImg);
    templateImg = single(templateImg);
    jointImg = single(jointImg);
    img = cat(4, subjImg, templateImg, jointImg);

    % Show detail viewer
    ea_imshowpair(img, options, title, callingfunction);
catch ME
    warning(ME.message);
    fprintf('Skip showing normalization...\n');
end

disp('Done.');
