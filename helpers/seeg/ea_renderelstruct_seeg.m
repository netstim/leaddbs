function [el_render, el_label, elSide] = ea_renderelstruct_seeg(options, resultfig, elstruct, pt, varargin)
% SEEG-specific rendering wrapper.
% Normal DBS rendering remains in the original ea_renderelstruct.m.

if ~exist('pt','var') || isempty(pt)
    pt = 1;
end

popts = options;

% Group-mode setup, matching the normal renderer where needed
if strcmp(options.leadprod,'group')
    [popts.root, popts.patientname] = fileparts(options.patient_list{pt});
    popts.root = [popts.root, filesep];

    recon = ea_regexpdir( ...
        [options.patient_list{pt}, filesep, 'reconstruction'], ...
        ['^', popts.patientname, '_desc-reconstruction\.mat$'], ...
        0, 'file');

    popts.subj.recon.recon = recon{1};
    popts = ea_detsides(popts);
end

elSide = popts.sides;

set(0, 'CurrentFigure', resultfig);

% SEEG renderer handles all electrode-coordinate cells itself
ea_trajectory_seeg(elstruct(pt), popts);

% Keep outputs compatible with ea_elvis
el_render = struct( ...
    'elpatch', [], ...
    'ellabel', [], ...
    'eltype', []);

% Preserve labels from previously rendered patients.
if ~isempty(varargin) && numel(varargin) >= 2
    el_label = varargin{2};
else
    el_label = [];
end

coords_mm = elstruct(pt).coords_mm;

% Retrieve shaft names, if available.
shaftNames = {};

if isfield(popts, 'subj') && ...
        isfield(popts.subj, 'recon') && ...
        isfield(popts.subj.recon, 'recon') && ...
        isfile(popts.subj.recon.recon)

    d = load(popts.subj.recon.recon, 'reco');

    if isfield(d, 'reco') && ...
            isfield(d.reco, 'props') && ...
            isfield(d.reco.props, 'elname')

        shaftNames = {d.reco.props.elname};
    end
end
% Use the axes belonging to the 3D figure.
ax = get(resultfig, 'CurrentAxes');

for i = 1:numel(coords_mm)

    coords = coords_mm{i};

    if isempty(coords) || size(coords,2) ~= 3
        continue
    end

    % Position label at the last contact.
    xyz = coords(end,:);

    % Determine shaft name.
    if ~isempty(shaftNames) && numel(shaftNames) >= i
        if iscell(shaftNames)
            shaftName = shaftNames{i};
        else
            shaftName = shaftNames(i);
        end
    else
        shaftName = sprintf('Shaft %d', i);
    end

    % Create label, initially hidden.
    h = text(ax, ...
        xyz(1), xyz(2), xyz(3), ...
        char(string(shaftName)), ...
        'Color', 'w', ...
        'FontSize', 15, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'Clipping', 'off', ...
        'Visible', 'off', ...
        'HitTest', 'off');

    el_label = [el_label; h];

end



end