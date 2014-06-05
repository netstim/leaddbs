function options=defaultoptions(options)


if ~isfield(options,'root')
options.root=input('Enter root directory...','s');
end
    
if ~isfield(options,'normalize')
options.normalize=0;
end

if ~isfield(options,'endtolerance')
options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
end

if ~isfield(options,'sprungwert')
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
end

if ~isfield(options,'refinesteps')
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
end

if ~isfield(options,'tra_stdfactor')
options.tra_stdfactor=0.9; % Default: 0.9 - the lower this factor, the lower the threshold (more included pixels in tra process).
end

if ~isfield(options,'cor_stdfactor')
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).
end

if ~isfield(options,'fit_cutoff')
options.fit_cutoff=1.8; % Cutoff for robust mean in trajectory reconstruction in standard deviations. Default 1.8.
end

if ~isfield(options,'verbose')
options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
end

if ~isfield(options,'sides')
options.sides=[1:2]; %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
end

if ~isfield(options,'maskwindow')
options.maskwindow=10; % size of the window that follows the trajectory
end

if ~isfield(options,'slow')
options.slow=0; % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.
end

if ~isfield(options,'axiscontrast')
options.axiscontrast=9; % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
end

if ~isfield(options,'zheights')
options.zheights=9; % if 1: use cor only, 2: use smoothed version of cor only, if 3: use mean of cor and tra, if 4: use multiplied cor * tra 5: use ^10 version of cor, 6: use ^10 version of cor.*tra 7: smoothed and then like 5. -> to determine heights of electrode contacts.
end

if ~isfield(options,'zresolution')
options.zresolution=10; % voxels are being parcellated into this amount of portions.
end

if ~isfield(options,'targetknown')
options.targetknown=0; % if 0, prior of stn will be used. Set priorstrength to 0 to not use any prior at all.
end

if ~isfield(options,'target')
options.target='stn';
end

if ~isfield(options,'priorstrength')
options.priorstrength=0.1; % how strong the prior information will be weighted. Set to 0 to not use prior at all.
end

if ~isfield(options,'showatlases')
options.showatlases=1;
end

if ~isfield(options,'showfibres')
options.showfibres=1;
end


if ~isfield(options,'showconnectivities')
options.showfibres=1;
end


if ~isfield(options,'writeoutimages')
options.writeoutimages=0;
end

if ~isfield(options,'overlays')
options.overlays={'edge',   'gp_h',     'gpi',  'morel'};
end

if ~isfield(options,'ocolors')
options.ocolors={'white',   'grey',    'red',   'multi'};       % use off or multi, red, blue, green, yellow, magenta, cyan, white, grey and black only.
end

if ~isfield(options,'oopacities')
options.oopacities=[1,    .2,         1,  	1];                 % use range 0-1
end

if ~isfield(options,'manualheightcorrection')
options.manualheightcorrection=1;
end

if ~isfield(options,'eldist')
options.eldist=2; % distance of electrode contacts. 3 mm for large (VIM) electrodes, 2 mm for small (STN/GPi) electrodes.
end