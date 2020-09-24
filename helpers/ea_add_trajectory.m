function ea_add_trajectory(~,~,options,pobj)


d=load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
if isfield(d.reco,'electrode')
    addedside=length(d.reco.electrode)+1;
else
    options=ea_detsides(options);
    addedside=max(options.sides)+1;
end

pobj.side=addedside;
pobj.hasMacro=0;
pobj.hasPlanning=1;
if exist('pobj','var')
    ea_trajectory(pobj);
end
