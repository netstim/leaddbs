function fv=ea_electrode2ply(directory,side,handles)

options=ea_handles2options(handles);
options.native=0; % MNI space only

if exist([directory,'ea_reconstruction.mat'],'file')
    load([directory,'ea_reconstruction.mat'],'reco');
    options.elmodel=ea_get_first_notempty_elmodel(reco.props);
end
options=ea_resolve_elspec(options);
[options.root,options.patientname]=fileparts(directory);
options.root=[options.root,filesep];
options.leadprod='dbs';
options.sidecolor=1;
options.prefs=ea_prefs;

[coords_mm,trajectory,markers]=ea_load_reconstruction(directory);
elstruct(1).coords_mm=coords_mm;
elstruct(1).trajectory=trajectory;
elstruct(1).name=options.patientname;
elstruct(1).markers=markers;
resultfig=figure('visible','off');

pobj.options=options;
pobj.elstruct=elstruct(1);
pobj.showMacro=1;
pobj.side=side;
set(0,'CurrentFigure',resultfig);
el_render(1)=ea_trajectory(pobj);

elrend=el_render.elpatch;

switch side
    case 1
        sidec='right_';
    case 2
        sidec='left_';
    otherwise
        sidec = num2str(side);
end

for f=1:length(elrend)
    fv(f).vertices=get(elrend(f),'Vertices');
    fv(f).faces=get(elrend(f),'Faces');
    fv(f).facevertexcdata=get(elrend(f),'FaceVertexCData');
end

fv=ea_concatfv(fv);
fv=ea_mapcolvert2ind(fv);
fv.faces=[fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)];

try
    plywrite([directory,'export',filesep,'ply',filesep,sidec,'electrode.ply'],fv.faces,fv.vertices,fv.facevertexcdata)
catch
    keyboard
end
