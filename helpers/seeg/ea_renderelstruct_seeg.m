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

el_label = [];

end