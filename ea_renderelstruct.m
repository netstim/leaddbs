function [el_render, el_label, elSide] = ea_renderelstruct(options,resultfig,elstruct,pt,el_render,el_label)
% Wrapper function to render lead trajectories based on elstruct

if ~exist('elstruct','var')
    [coords_mm,trajectory,markers] = ea_load_reconstruction(options);
    elstruct(1).coords_mm = coords_mm;
    elstruct(1).trajectory = trajectory;
    elstruct(1).name = options.patientname;
    elstruct(1).markers = markers;
end

if ~exist('pt','var')
    pt = 1;
end

popts = options;
if strcmp(options.leadprod,'group')
    [popts.root, popts.patientname] = fileparts(options.patient_list{pt});
    popts.root = [popts.root, filesep];
    recon = ea_regexpdir([options.patient_list{pt}, filesep, 'reconstruction'], ['^', popts.patientname,'_desc-reconstruction\.mat$'], 0, 'file');
    popts.subj.recon.recon = recon{1};
    popts = ea_detsides(popts);
end

elSide = popts.sides;

for side=elSide
    try
        pobj = ea_load_electrode(options.subj.recon.recon, side);
        pobj.hasPlanning = 1;
        pobj.showPlanning = strcmp(options.leadprod,'or');
    end
    pobj.pt = pt;
    pobj.options = popts;
    pobj.elstruct = elstruct(pt);
    pobj.showMacro = 1;
    pobj.side = side;

    set(0,'CurrentFigure',resultfig);
    if ~exist('el_render','var') || ~isstruct(el_render)
        el_render = struct([]);
    end
    options.leadprod = 'seeg';
    if isequal(options.leadprod, 'seeg')
        ea_trajectory_seeg(elstruct, options);
        stub = struct('ellabel', 'SEEG');
        if isempty(el_render)
            el_render = stub;
        else
            el_render(end+1) = stub;
        end
    else
        if isempty(el_render)
            el_render = ea_trajectory(pobj);
        else
            el_render(end+1) = ea_trajectory(pobj);
        end
    end
    
    if ~exist('el_label','var')
        if isfield(el_render(end),'ellabel') && ~isempty(el_render(end).ellabel)
            el_label = el_render(end).ellabel;
        else
            el_label = 'SEEG';
        end
    else
        try
            el_label(end+1) = el_render(end).ellabel;
        catch
            if iscell(el_label)
                el_label{end+1} = 'SEEG';
            else
                el_label = [el_label, 'SEEG'];
            end
        end
    end
end
