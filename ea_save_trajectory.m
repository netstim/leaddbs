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
    try trajectories(obj.site).dbs.reco.(spacenames{spacenm}).coords_mm{obj.site}=reco.(spacenames{spacenm}).coords_mm{obj.site}; end
    try trajectories(obj.site).dbs.reco.(spacenames{spacenm}).trajectory{obj.site}=reco.(spacenames{spacenm}).trajectory{obj.site}; end
    try trajectories(obj.site).dbs.reco.(spacenames{spacenm}).markers(obj.site)=reco.(spacenames{spacenm}).markers(obj.site); end
end
trajectories(obj.site).dbs.elmodel=obj.elmodel;
%%%

trajectories(obj.site).planning.color=obj.color;
trajectories(obj.site).planning.planRelative=obj.planRelative;
trajectories(obj.site).planning.target=obj.target;

trajectories(obj.site).micro.relateMicro=obj.relateMicro;
save([obj.options.root,obj.options.patientname,filesep,'ea_trajectories.mat'],'trajectories');

