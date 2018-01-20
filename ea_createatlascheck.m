function ea_createatlascheck(h,gf,options)
%     options.d2=options.prefs.machine.d2;
options.atlasset=h.Parent.Label;
load([ea_space(options,'atlases'),options.atlasset,filesep,'atlas_index.mat']);
[~,six]=ismember(h.Label,ea_rmext(atlases.names));
fv=atlases.fv(six,:);
pixdim=atlases.pixdim(six,:);

if length(fv)>1
   xyz=[fv{1}.vertices;fv{2}.vertices];
   pixdim=mean([pixdim{1};pixdim{2}]);
else
    xyz=fv{1}.vertices;
   pixdim=pixdim{1};
end
mz=mean(xyz);
vmz=abs(max(xyz)-min(xyz));
options.d2.writeatlases=1;
options.d2.col_overlay=0;
options.d2.con_color=[0.8,0.1,0];
options.d2.atlasopacity=0.2;
options.d2.tracor=1;
options.d2.bbsize=(max(vmz)/1.7);
options.d2.depth=mz;
options.d2.showlegend=0;
options.d2.showstructures={h.Label};
[Vtra,Vcor,Vsag]=ea_assignbackdrop(['Patient Pre-OP (',h.Parent.Parent.Label,')'],options,'Patient');
Vs={Vtra,Vcor,Vsag};
options.sides=1;
evalin('base','custom_cont=2;');
hf=ea_writeplanes(options,options.d2.depth,options.d2.tracor,Vs{options.d2.tracor},'on',2);
set(hf,'Position',[0,0,800,800]);