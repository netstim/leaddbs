function ea_save_electrode(obj)

if isfile(obj.options.subj.recon.recon)
    d = load(obj.options.subj.recon.recon);
else
    d.reco.electrode = struct;
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

d.reco.electrode(obj.side).micro.relateMicro=obj.relateMicro;
save(obj.options.subj.recon.recon, '-struct', 'd');

