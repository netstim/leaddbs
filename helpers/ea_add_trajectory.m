function ea_add_trajectory(~,~,options,pobj,side)

d=load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
if ~exist('side', 'var')
    if isfield(d.reco,'electrode')
        side=length(d.reco.electrode)+1;
    else
        options=ea_detsides(options);
        side=max(options.sides)+1;
    end
end

pobj.side = side;
pobj.hasMacro=0;
pobj.hasPlanning=1;

try
    pobj.planningAppearance=d.reco.electrode(pobj.side).plan.planningAppearance;
end
try
    pobj.plan2elstruct=d.reco.electrode(pobj.side).plan.plan2elstruct;
end
try
    pobj.plan2elstruct_model=d.reco.electrode(pobj.side).plan.plan2elstruct_model;
end

if exist('pobj','var')
    ea_trajectory(pobj);
end
