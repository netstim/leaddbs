function ea_show_coregistration(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

disp('Preparing images to show Coregistration...');

% Load & adjust coregistered and anchor images
moving = ea_load_nii(options.moving);
fixed = ea_load_nii([options.fixed]);

movingImg = moving.img;
movingImg(:) = ea_nanzscore(movingImg(:));
movingImg = (movingImg+2.5)/5;
movingImg(movingImg<0) = 0;
movingImg(movingImg>1) = 1;

fixedImg = fixed.img;
fixedImg(:) = ea_nanzscore(fixedImg(:));
fixedImg = (fixedImg+2.5)/5;
fixedImg(fixedImg<0) = 0;
fixedImg(fixedImg>1) = 1;

% Set joint img
jointImg = cat(4, 0.1*fixedImg+0.9*movingImg, 0.4*fixedImg+0.6*movingImg, 0.9*fixedImg+0.1*movingImg);

% Set image for detail viewer
img = cat(4, fixedImg, movingImg, jointImg);

% Show detail viewer
ea_imshowpair(img,options,options.tag,'coregistration');

disp('Done.');
