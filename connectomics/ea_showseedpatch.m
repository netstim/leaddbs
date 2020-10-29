function [seedsurf]=ea_showseedpatch(resultfig,pV,sX,options)

set(0,'CurrentFigure',resultfig)
bb=[1,1,1;size(sX)];
bb=ea_vox2mm(bb,pV.mat);
gv=cell(3,1);
for dim=1:3
    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(sX,dim));
end
[X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});


fvs=isosurface(X,Y,Z,permute(sX,[2,1,3]),0.5); % seed

if ischar(options.prefs.lhullsimplify)

    simplify=1000/length(fvs.faces);
    fvs=reducepatch(fvs,simplify);

else
    if options.prefs.lhullsimplify<1 && options.prefs.lhullsimplify>0
        fvs=reducepatch(fvs,options.prefs.lhullsimplify);
    end
end

seedsurf=patch(fvs,'FaceColor',options.prefs.lc.seedsurfc,'facealpha',0.7,'EdgeColor','none','facelighting','phong');
set(seedsurf,'DiffuseStrength',0.9)
set(seedsurf,'SpecularStrength',0.1)
set(seedsurf,'FaceAlpha',0.3);
