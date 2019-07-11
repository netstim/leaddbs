function pobj=ea_load_trajectory(directory,side)
% 
%
% USAGE:
%
%    pobj = ea_load_trajectory(directory,side)
%
% INPUTS:
%    directory:
%    side:
%
% OUTPUT:
%    pobj:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ningfei Li, Original file
%       - Daniel Duarte, Documentation

load([directory,'ea_trajectories.mat']);

pobj.elmodel=trajectories(side).dbs.elmodel;
%%%
pobj.color=trajectories(side).planning.color;
pobj.planRelative=trajectories(side).planning.planRelative;
pobj.target=trajectories(side).planning.target;
pobj.relateMicro=trajectories(side).micro.relateMicro;

