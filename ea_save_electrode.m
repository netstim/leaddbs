function ea_save_electrode(obj)

if exist([obj.options.root,obj.options.patientname,filesep,'ea_reconstruction.mat'],'file')
    d=load([obj.options.root,obj.options.patientname,filesep,'ea_reconstruction.mat']);
else
    d.reco.electrode=struct;
end


spacenames={'native','scrf','mni','acpc'};
for spacenm=1:length(spacenames)
    try d.reco.electrode(obj.side).dbs.reco.(spacenames{spacenm}).coords_mm{obj.side}=reco.(spacenames{spacenm}).coords_mm{obj.side}; end
    try d.reco.electrode(obj.side).dbs.reco.(spacenames{spacenm}).trajectory{obj.side}=reco.(spacenames{spacenm}).trajectory{obj.side}; end
    try d.reco.electrode(obj.side).dbs.reco.(spacenames{spacenm}).markers(obj.side)=reco.(spacenames{spacenm}).markers(obj.side); end
end
d.reco.electrode(obj.side).dbs.elmodel=obj.elmodel;
%%%

d.reco.electrode(obj.side).plan.color=obj.color;
d.reco.electrode(obj.side).plan.planRelative=obj.planRelative;
d.reco.electrode(obj.side).plan.target=obj.target;
d.reco.electrode(obj.side).plan.planningAppearance=obj.planningAppearance;
d.reco.electrode(obj.side).plan.plan2elstruct=obj.plan2elstruct;
d.reco.electrode(obj.side).plan.plan2elstruct_model=obj.plan2elstruct_model;

d.reco.electrode(obj.side).micro.relateMicro=obj.relateMicro;
save([obj.options.root,obj.options.patientname,filesep,'ea_reconstruction.mat'],'-struct','d');

