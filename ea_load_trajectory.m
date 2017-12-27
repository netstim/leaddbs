function pobj=ea_load_trajectory(directory,side)

    load([directory,'ea_trajectories.mat']);

pobj.elmodel=trajectories(side).dbs.elmodel;
%%%
pobj.color=trajectories(side).planning.color;
pobj.planRelative=trajectories(side).planning.planRelative;
pobj.target=trajectories(side).planning.target;
pobj.relateMicro=trajectories(side).micro.relateMicro;

