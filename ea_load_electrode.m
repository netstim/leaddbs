function pobj=ea_load_electrode(directory,side)

d=load([directory,'ea_reconstruction.mat']);

pobj.elmodel=d.reco.electrode(side).dbs.elmodel;
%%%
pobj.color=d.reco.electrode(side).plan.color;
pobj.planRelative=d.reco.electrode(side).plan.planRelative;
pobj.target=d.reco.electrode(side).plan.target;
pobj.relateMicro=d.reco.electrode(side).micro.relateMicro;

