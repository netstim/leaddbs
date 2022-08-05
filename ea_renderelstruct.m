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
    if exist('el_render','var')
        el_render(end+1) = ea_trajectory(pobj);
    else
        el_render(1) = ea_trajectory(pobj);
    end

    if ~exist('el_label','var')
        el_label = el_render(end).ellabel;
    else
        try
            el_label(end+1) = el_render(end).ellabel;
        end
    end
end
