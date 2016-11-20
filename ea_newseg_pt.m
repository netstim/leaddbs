function ea_newseg_pt(options,dartel)

% wrapper for ea_newseg multimodal to be called whenever generating
% segmentations on preop MRI

directory=[options.root,options.patientname,filesep];
ea_checkcoregallmri(options,0,0); % dont use FA here. Except that may borrow same function as in ANTs multimodal.
[options,presentfiles]=ea_assignpretra(options);
ea_newseg(directory,presentfiles,dartel,options,1)