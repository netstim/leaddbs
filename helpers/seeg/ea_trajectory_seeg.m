function ea_trajectory_seeg(elstruct, options)
% Render SEEG contacts and curved electrode models.
%
% Contact spheres use the original LeGUI rendering style.
% The full electrode shaft and contacts are rendered using
% ea_render_curved_lead.

if ~isfield(elstruct, 'coords_mm')
    error('Input must have field coords_mm (1×N cell of [n×3]).');
end
% Choose how SEEG electrodes are displayed.
renderMode = 'Electrode model';

% Read the selection saved by the 2D/3D Settings window.
if isfield(options, 'prefs') && ...
        isfield(options.prefs, 'machine') && ...
        isfield(options.prefs.machine, 'd2') && ...
        isfield(options.prefs.machine.d2, 'seegRenderMode') && ...
        ~isempty(options.prefs.machine.d2.seegRenderMode)

    renderMode = ...
        options.prefs.machine.d2.seegRenderMode;
end

showSpheres = strcmpi(renderMode, 'Spheres');

showElectrodeModel = ...
    strcmpi(renderMode, 'Electrode model');
% Sphere resolution and size
[Xs, Ys, Zs] = sphere(12);
r_sphere = 1.2;  % mm

hold on

Ncells = numel(elstruct.coords_mm);
rendered_any = false;

% Same shaft color palette used by LeGUI Assign/Check mode
nShafts = Ncells;

hues = linspace(0, 1, nShafts + 1)';
hues(end) = [];

hues = mod(hues + 0.18, 1);

shaftColors = hsv2rgb([ ...
    hues, ...
    0.65*ones(nShafts,1), ...
    0.85*ones(nShafts,1)]);

% Avoid red-like colors because red is typically used for highlighting
redLike = shaftColors(:,1) > 0.75 & ...
    shaftColors(:,2) < 0.45 & ...
    shaftColors(:,3) < 0.45;

shaftColors(redLike,:) = ...
    repmat([0.00 0.60 0.80], sum(redLike), 1);
for i = 1:Ncells

    coords = elstruct.coords_mm{i};

    if isempty(coords) || size(coords,2) ~= 3
        continue
    end

    % Remove exact duplicate locations for purposes of determining whether
    % a shaft direction can be defined. Original contact coordinates are untouched.
    uniqueCoords = unique(coords, 'rows', 'stable');

    if showElectrodeModel && size(uniqueCoords,1) < 2
        warning(['SEEG shaft %d has only one unique contact position. ' ...
            'Cannot determine electrode direction; skipping model.'], i);
        continue
    end
    if showSpheres && ...
            isfield(options, 'subjectColor') && ...
            ~isempty(options.subjectColor)

        % All shafts/contacts from this subject use the same color
        ccol = options.subjectColor;

    else

        % Original per-shaft coloring
        ccol = shaftColors(i,:);

    end

    if showSpheres

        % Draw the original LeGUI contact spheres.
        for k = 1:size(coords,1)

            x = coords(k,1);
            y = coords(k,2);
            z = coords(k,3);

            surf( ...
                r_sphere*Xs + x, ...
                r_sphere*Ys + y, ...
                r_sphere*Zs + z, ...
                'FaceColor', ccol, ...
                'EdgeColor', 'none');
        end
    end

    % Build the curved centerline.
    if isfield(elstruct, 'curve_mm') && ...
            numel(elstruct.curve_mm) >= i && ...
            ~isempty(elstruct.curve_mm{i}) && ...
            size(elstruct.curve_mm{i},2) == 3

        % Use an explicitly saved curve when available.
        Ct = elstruct.curve_mm{i};

    elseif size(coords,1) >= 3

        % Fit a quadratic curve through the first, middle,
        % and last contact.
        P0 = coords(1,:);
        P1 = coords(end,:);
        Pm = coords(round((size(coords,1)+1)/2),:);

        A = [0.25, 0.5; 1, 1];
        rhs = [Pm-P0; P1-P0];

        xAB = A \ rhs;

        a = xAB(1,:);
        b = xAB(2,:);
        c0 = P0;

        t = linspace(0,1,200).';

        Ct = a.*(t.^2) + b.*t + c0;

    else

        % Straight-line fallback.
        t = linspace(0,1,100).';

        Ct = ...
            (1-t).*coords(1,:) + ...
            t.*coords(end,:);

    end

    if showElectrodeModel

        % Resolve the selected electrode model and render settings.
        [elspec, renderopt] = ...
            ea_resolve_curved_seeg_rendering( ...
            elstruct, options, i);

        % Render the actual curved SEEG electrode.
        ea_render_curved_lead( ...
            Ct, ...
            coords, ...
            elspec, ...
            renderopt);
    end
    rendered_any = true;
end

hold off

if ~rendered_any
    return
end
end


function [elspec, renderopt] = ...
    ea_resolve_curved_seeg_rendering( ...
    elstruct, options, side)

elspec = struct;
renderopt = struct;

elmodel = '';

% First choice: model stored in elstruct.props.
if isfield(elstruct, 'props') && ...
        numel(elstruct.props) >= side && ...
        isfield(elstruct.props, 'elmodel') && ...
        ~isempty(elstruct.props(side).elmodel)

    elmodel = elstruct.props(side).elmodel;

    % Second choice: model stored directly in elstruct.
elseif isfield(elstruct, 'elmodel') && ...
        ~isempty(elstruct.elmodel)

    if iscell(elstruct.elmodel) && ...
            numel(elstruct.elmodel) >= side

        elmodel = elstruct.elmodel{side};

    elseif isstring(elstruct.elmodel) || ...
            ischar(elstruct.elmodel)

        elmodel = char(elstruct.elmodel);
    end

    % Final fallback: model in options.
elseif isfield(options, 'elmodel') && ...
        ~isempty(options.elmodel)

    if iscell(options.elmodel) && ...
            numel(options.elmodel) >= side

        elmodel = options.elmodel{side};

    else
        elmodel = options.elmodel;
    end
end

if isstring(elmodel)
    elmodel = char(elmodel);

elseif iscell(elmodel) && ~isempty(elmodel)
    elmodel = elmodel{1};
end

% Resolve model geometry.
if ~isempty(elmodel) && ...
        exist('ea_resolve_elspec', 'file') == 2

    try
        specopts = struct('elmodel', elmodel);
        specopts = ea_resolve_elspec(specopts);

        if isfield(specopts, 'elspec')
            elspec = specopts.elspec;
        end
    catch
    end
end

% Fall back to options.elspec.
if isempty(fieldnames(elspec)) && ...
        isfield(options, 'elspec')

    elspec = options.elspec;
end

% Active contacts.
if isfield(elstruct, 'activecontacts') && ...
        numel(elstruct.activecontacts) >= side && ...
        ~isempty(elstruct.activecontacts{side})

    renderopt.active_idx = ...
        find(elstruct.activecontacts{side});

else
    renderopt.active_idx = [];
end

if isfield(options, 'd3') && ...
        isfield(options.d3, 'hlactivecontacts') && ...
        ~options.d3.hlactivecontacts

    renderopt.active_idx = [];
end

% Transparency settings.
if isfield(options, 'd3') && ...
        isfield(options.d3, 'elrendering') && ...
        options.d3.elrendering == 2

    renderopt.contact_alpha = 0.1;
    elspec.insulation_alpha = 0.1;

else
    renderopt.contact_alpha = 1;
    elspec.insulation_alpha = 1;
end

renderopt.faceLighting = 'phong';
renderopt.specularStrength = 0.2;
renderopt.contactSpecularStrength = 1.0;
end