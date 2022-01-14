function ea_imoverlay(refNifti, foregroundNifti, type)
% Overlay two NIfTI images as in the "Check Registration" window.

if ~exist('type', 'var') || isempty(type)
    type = 'coregistration';
else
    switch type
        case {'coreg', 'coregistration'}
            type = 'coregistration';
        case {'norm', 'normalization', 'normalization dbs'}
            type = 'normalization dbs';
    end
end

ref = ea_load_nii(refNifti);
foreground = ea_load_nii(foregroundNifti);

[~, foregroundName] = ea_niifileparts(foreground.fname);
[~, refName] = ea_niifileparts(ref.fname);

tag = ['Foreground: ', foregroundName, ' & Reference: ', refName];
options.patientname = '';

switch type
    case 'coregistration'
        foreground.img(:) = ea_nanzscore(foreground.img(:));
        foreground.img = (foreground.img + 2.5)/5; % set max/min to +/- 2.5 standard deviations
        foreground.img(foreground.img<0) = 0;
        foreground.img(foreground.img>1) = 1;

        ref.img(:) = ea_nanzscore(ref.img(:));
        ref.img = (ref.img + 2.5)/5; % set max/min to +/- 2.5 standard deviations
        ref.img(ref.img<0) = 0;
        ref.img(ref.img>1) = 1;

        jointImg = cat(4, 0.1*ref.img+0.9*foreground.img, 0.4*ref.img+0.6*foreground.img, 0.9*ref.img+0.1*foreground.img);

        Images = cat(4, foreground.img, ref.img, jointImg);

        helptext = {'Click to show reference image', ...
                    'Use </> to decrease/increase box size while clicking', ...
                    'Arrow keys / Mouse wheel: Scroll through image', '', ...
                    'Z: Zoom in/out', ...
                    'X: Hybrid view on/off', ...
                    'A: Axial view', ...
                    'C: Coronal view', ...
                    'S: Saggital view'};
        ea_imshowpair(Images, options, tag, 'coregistration', helptext);
    case 'normalization dbs'
        foreground.img = (foreground.img - min(foreground.img(:))) / max(foreground.img(:));
        foreground.img(foreground.img>0.5) = 0.5;
        foreground.img = (foreground.img - min(foreground.img(:))) / max(foreground.img(:));

        ref.img(:) = zscore(ref.img(:));
        ref.img = (ref.img-min(ref.img(:)))/(max(ref.img(:))-min(ref.img(:)));

        load([ea_space, 'wires.mat'], 'wires');
        wires = single(wires);
        wires = wires/255;
        wires = wires*0.2;
        wires = wires+0.8;

        jointImg = foreground.img .* wires;
        jointImgStd = std(jointImg(:));
        jointImg = jointImg - min(jointImg(:));
        jointImg = jointImg ./ (max(jointImg(:)) - jointImgStd);
        jointImg(jointImg>1) = 1;

        foreground.img = single(foreground.img);
        ref.img = single(ref.img);
        jointImg = single(jointImg);

        Images = cat(4, foreground.img, ref.img, jointImg);

        helptext = {'Click to show reference image', ...
                    'Use </> to decrease/increase box size while clicking', ...
                    'Arrow keys / Mouse wheel: Scroll through image','', ...
                    'Alt+1,2,...: Switch between available templates [FA=0]', ...
                    'Z: Zoom in/out', ...
                    'X: Hybrid view on/off', ...
                    'A: Axial view', ...
                    'C: Coronal view', ...
                    'S: Saggital view'};
        ea_imshowpair(Images, options, tag, 'normalization dbs', helptext);
end
