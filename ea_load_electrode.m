function pobj=ea_load_electrode(recon_path,side)

d=load(recon_path);

pobj.elmodel=d.reco.electrode(side).dbs.elmodel;
pobj.color=d.reco.electrode(side).plan.color;
pobj.planRelative=d.reco.electrode(side).plan.planRelative;
pobj.target=d.reco.electrode(side).plan.target;
pobj.relateMicro=d.reco.electrode(side).micro.relateMicro;
