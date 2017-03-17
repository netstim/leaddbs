function fv=ea_electrode2ply(directory,side,handles)

options=ea_handles2options(handles);
options=ea_resolve_elspec(options);
[options.root,options.patientname]=fileparts(directory);
options.root=[options.root,filesep];



[coords_mm,trajectory,markers]=ea_load_reconstruction(options);

elstruct(1).coords_mm=coords_mm;
elstruct(1).coords_mm=ea_resolvecoords(markers,options);
elstruct(1).trajectory=trajectory;
elstruct(1).name=options.patientname;
elstruct(1).markers=markers;
resultfig=figure('visible','off');

[el_render(1).el_render,el_label(:,1)]=ea_showelectrode(resultfig,elstruct(1),1,options);
elrend=el_render.el_render;

    switch side
        case 1
            sidec='right_';
        case 2
            sidec='left_';
    end
    for f=1:length(elrend{side})
        fv(f).vertices=get(elrend{side}(f),'Vertices');
        fv(f).faces=get(elrend{side}(f),'Faces');
        fv(f).facevertexcdata=get(elrend{side}(f),'FaceVertexCData');
    end
    fv=ea_concatfv(fv);
    fv=ea_mapcolvert2ind(fv);
fv.faces=[fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)];

ea_patch2ply([directory,'export',filesep,'ply',filesep,sidec,'electrode.ply'],fv.vertices',fv.faces',fv.facevertexcdata');



