function ea_save_trajectory(obj)

if exist([obj.options.root,obj.options.patientname,filesep,'ea_trajectories.mat'],'file');
    load([obj.options.root,obj.options.patientname,filesep,'ea_trajectories.mat']);
else
    trajectories=struct;
end
if exist([obj.options.root,obj.options.patientname,filesep,'ea_reconstruction.mat'],'var')
    load([obj.options.root,obj.options.patientname,filesep,'ea_reconstruction.mat']);
end


spacenames={'native','scrf','mni','acpc'};
for spacenm=1:length(spacenames)
    try trajectories(obj.side).dbs.reco.(spacenames{spacenm}).coords_mm{obj.side}=reco.(spacenames{spacenm}).coords_mm{obj.side}; end
    try trajectories(obj.side).dbs.reco.(spacenames{spacenm}).trajectory{obj.side}=reco.(spacenames{spacenm}).trajectory{obj.side}; end
    try trajectories(obj.side).dbs.reco.(spacenames{spacenm}).markers(obj.side)=reco.(spacenames{spacenm}).markers(obj.side); end
end
trajectories(obj.side).dbs.elmodel=obj.elmodel;
%%%

trajectories(obj.side).planning.color=obj.color;
trajectories(obj.side).planning.planRelative=obj.planRelative;
trajectories(obj.side).planning.target=obj.target;

trajectories(obj.side).micro.relateMicro=obj.relateMicro;
save([obj.options.root,obj.options.patientname,filesep,'ea_trajectories.mat'],'trajectories');

