function pobj=ea_load_trajectory(directory,site)

    load([directory,'ea_trajectories.mat']);

pobj.elmodel=trajectories(site).dbs.elmodel;
%%%
pobj.color=trajectories(site).planning.color;
pobj.planRelative=trajectories(site).planning.planRelative;
pobj.target=trajectories(site).planning.target;
pobj.relateMicro=trajectories(site).micro.relateMicro;

