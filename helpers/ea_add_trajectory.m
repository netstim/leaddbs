function ea_add_trajectory(~,~,options,pobj,side)

d=load([options.root,options.patientname,filesep,'reconstruction',filesep,options.patientname,'_desc-reconstruction.mat']);
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
if exist('pobj','var')
    ea_trajectory(pobj);
end
