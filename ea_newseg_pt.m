function ea_newseg_pt(options,dartel,del,force)

% wrapper for ea_newseg multimodal to be called whenever generating
% segmentations on preop MRI
if ~exist('force','var')
    force=0;
end

directory=[options.root,options.patientname,filesep];
ea_checkcoregallmri(options,0,0); % dont use FA here. Except that may borrow same function as in ANTs multimodal.
[options,presentfiles]=ea_assignpretra(options);
ea_newseg(directory,presentfiles,dartel,options,del,force)