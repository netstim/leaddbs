function options=ea_defaultoptions(varargin)
try 
    options=varargin{1};
catch
    options=struct;
end

if ~isfield(options,'normalize')
    options.normalize=0;
end

if ~isfield(options,'endtolerance')
    options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
end

if ~isfield(options,'sprungwert')
    options.sprungwert=4; % how far electrode centroid may lie (in xy axis) from last to current slice.
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
    options.sides=[1:2]; %side=1 -> right electrode, side=2 -> left electrode. both: [1:2]
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


if ~isfield(options,'d2')
    options.d2.write=1;
    options.d2.atlasopacity=0.5;
    options.d2.col_overlay=1;
    options.d2.con_overlay=1;
    options.d2.lab_overlay=1;
    options.d2.bbsize=50;
end

if ~isfield(options.d2,'write')
    options.d2.write=1;
end

if ~isfield(options.d2,'atlasopacity')
    options.d2.atlasopacity=0.5;
end

if ~isfield(options.d2,'col_overlay')
    options.d2.col_overlay=1;
end

if ~isfield(options.d2,'con_overlay')
    options.d2.con_overlay=1;
end

if ~isfield(options.d2,'lab_overlay')
    options.d2.lab_overlay=1;
end

if ~isfield(options.d2,'bbsize')
    options.d2.bbsize=100;
end


if ~isfield(options,'d3')
    options.d3.write=1;
    options.d3.prolong_electrode=2;
    options.d3.writeatlases=1;
end

if ~isfield(options.d3,'write')
    options.d2.write=1;
end


if ~isfield(options.d3,'prolong_electrode')
    options.d2.prolong_electrode=1;
end

if ~isfield(options.d3,'writeatlases')
    options.d2.writeatlases=1;
end

if ~isfield(options.d3,'hlactivecontacts')
    options.d3.hlactivecontacts=0;
end

if ~isfield(options.d3,'elrendering')
    options.d3.elrendering=1;
end

if ~isfield(options,'dostimulation')
    options.dostimulation=1;
end

if ~isfield(options,'writeoutimages')
    options.writeoutimages=0;
end

if ~isfield(options,'refinelocalization')
    options.refinelocalization=1;
end

if ~isfield(options,'elmodel')
    options.elmodel='Medtronic 3389'; % Specify electrode model here. Review or add available electrode models in ea_resolve_elspec.m
end

if ~isfield(options,'earoot')
    
    options.earoot=ea_getearoot; 
end

if ~isfield(options,'fiberthresh')
    
    options.fiberthresh=5; % 
end

if ~isfield(options,'writeoutstats')
    
    options.writeoutstats=1;
end

if ~isfield(options,'labelatlas')
    
    options.labelatlas='aal';
end


if ~isfield(options,'atlasset')
    
    options.atlasset='STN_GPi';
end

if ~isfield(options,'writeoutpm')
    
    options.writeoutpm = 0;
end

if ~isfield(options,'prefs')
    options.prefs=ea_prefs;
end
