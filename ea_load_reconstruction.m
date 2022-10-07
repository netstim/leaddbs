function  [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(varargin)
% [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(varargin)

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
%
% Andreas Horn
%
% Modified for groupmode 04/2017 by Ari Kappel
% Manually reversed uipatdirs changes 08/26/17 Andy Horn
% Please do not use uipatdirs to determine patient directory, this will
% confuse calls with multiple patients selected.

if isstruct(varargin{1})
    options = varargin{1};
else % TODO: check where this function is called and see how inputs are formated
    directory=varargin{1};
    if ~strcmp(directory(end),filesep)
        directory=[directory,filesep];
    end
    options=ea_getptopts(directory);
end

if ~isfield(options, 'subj')
    options.subj.recon.recon = [options.root,options.patientname,filesep,'reconstruction',filesep,options.patientname,'_desc-reconstruction.mat'];
end

try
    % Load Reconstruction
    load(options.subj.recon.recon, 'reco');
catch
    ea_cprintf('CmdWinWarnings', 'Failed to load reconstruction for %s!\n', options.subj.subjId);
end

if ~isfield(reco,'native') && isfield(reco,'mni') && options.native
    ea_reconstruction2native(options);
    load(options.subj.recon.recon);
elseif isfield(reco,'native') && ~isfield(reco,'mni') && ~options.native
    ea_reconstruction2mni(options);
    load(options.subj.recon.recon);
end

if options.native
    if isfield(options, 'loadnativereco') && options.loadnativereco || ~isfield(reco, 'scrf')
        % Load native reco for recalculation or manual reconstruction
        space_type = 'native';
    else
        % Load scrf reco by default
        space_type = 'scrf';
    end
else
    space_type = 'mni';
end

markers = reco.(space_type).markers;
if ~isfield(markers,'x')
    for side=1:2
        [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
        markers(side).x = markers(side).head+xunitv*(options.elspec.lead_diameter/2);
        markers(side).y = markers(side).head+yunitv*(options.elspec.lead_diameter/2);
    end
end

if isfield(reco.(space_type),'coords_mm')
    coords_mm = reco.(space_type).coords_mm;
else
    [coords_mm]=ea_resolvecoords(markers,options,0);
end

if isfield(reco.(space_type),'trajectory')
    trajectory = reco.(space_type).trajectory;
else
    [~,trajectory,markers]=ea_resolvecoords(markers,options,0);
end

try
    coords_acpc=reco.acpc.coords_mm;
catch
    coords_acpc=nan;
end

try
    manually_corrected=reco.props(options.elside).manually_corrected;
    elmodel=reco.props(options.elside).elmodel;
    % Fix Abbott lead name
    elmodel = strrep(elmodel, 'St. Jude', 'Abbott');
catch % legacy
    [elmodel,first_notempty_side]=ea_get_first_notempty_elmodel(reco.props);
    manually_corrected=reco.props(first_notempty_side).manually_corrected;
end

%if elmodel is empty, search for the first available side that has a model
if isempty(elmodel)
    for side=1:length(reco.props)
        elmodel=reco.props(side).elmodel;
        if ~isempty(elmodel)
            break
        end
    end
end

if ~exist('manually_corrected','var')
    manually_corrected=0;
end
if isempty(manually_corrected)
    manually_corrected=0;
end

if ~exist('elmodel','var')
    elmodel=[];
end
