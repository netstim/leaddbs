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

coords_acpc = nan; % initialize output

if isstruct(varargin{1})
    options = varargin{1};
else % TODO: check where this function is called and see how inputs are formated
    directory=varargin{1};
    if ~strcmp(directory(end),filesep)
        directory=[directory,filesep];
    end
    options=ea_getptopts(directory);
end

try
    % Load Reconstruction
    load(options.subj.recon.recon, 'reco');
end

if exist('reco','var')
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

else % legacy format
    if ~isfield(options,'elspec')
        options=ea_getptopts(directory,options);
    end
    if ~exist('markers','var') % backward compatibility to old recon format
        for iside=1:length(options.sides)
            side=options.sides(iside);

            markers(side).head=coords_mm{side}(1,:);
            markers(side).tail=coords_mm{side}(4,:);
            [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
            markers(side).x = coords_mm{side}(1,:)+xunitv*(options.elspec.lead_diameter/2);
            markers(side).y = coords_mm{side}(1,:)+yunitv*(options.elspec.lead_diameter/2);
        end

        elmodel=options.elmodel;
        % this is still legacy format but has markers variable in it now.
        save([options.root,options.patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm','markers','elmodel');
    end

    try
        load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
    catch % generate trajectory from coordinates.
        trajectory{1}=ea_fit_line(coords_mm(1:4,:));
        trajectory{2}=ea_fit_line(coords_mm(options.elspec.numel+1:options.elspec.numel+4,:));
    end

    % we have all variables now but need to put them into reco format and
    % they are in MNI.
    options.native=0;
    options.hybridsave=1;

    try
        ea_save_reconstruction(coords_mm,trajectory,markers,'Medtronic 3389',0,options);
        options.native=1;
        [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
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
