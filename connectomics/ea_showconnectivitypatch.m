function [matsurf,labels]=ea_showconnectivitypatch(resultfig,pV,mX,thresh,atlaslegend,atlasindices,usecolormap,showregs,showlabels)
matsurf=[];
labels=[];

if showregs
    tmX=abs(mX)>thresh;

    options.prefs=ea_prefs('');
    if options.prefs.hullsmooth
        tmX=smooth3(tmX,'gaussian',options.prefs.hullsmooth);
    end

    bb=[1,1,1;size(mX)];
    bb=ea_vox2mm(bb,pV.mat);
    gv=cell(3,1);
    for dim=1:3
        gv{dim}=linspace(bb(1,dim),bb(2,dim),size(mX,dim));
    end
    [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
    fv=isosurface(X,Y,Z,permute(tmX,[2,1,3]),0.5); % graph_metric

    if ischar(options.prefs.hullsimplify)
        % get to 1500 faces
        simplify=1500;
        fv=reducepatch(fv,simplify);
    else
        if options.prefs.hullsimplify<1 && options.prefs.hullsimplify>0
            fv=reducepatch(fv,options.prefs.hullsimplify);
        elseif options.prefs.hullsimplify>1
            simplify=options.prefs.hullsimplify/length(fv.faces);
            fv=reducepatch(fv,simplify);
        end
    end
    try fv=ea_smoothpatch(fv,[],20); end

    set(0,'CurrentFigure',resultfig)
    matsurf=patch(fv,'facealpha',0.7,'EdgeColor','none','facelighting','phong','FaceColor','interp');

    zidx=mX==0;
    zidx=logical(zidx+isnan(mX));

    % %cgX(isnan(cgX))=0;
    % try % if cgX isn't empty..
    % cgX(cgX~=0)=cgX(cgX~=0)-ea_nanmin(cgX(cgX~=0));
    %
    % cgX(cgX~=0)=(cgX(cgX~=0)/ea_nanmax(cgX(cgX~=0)))*64;
    % end

    mX(zidx)=0; % reset prior zero/nan values to zero

    nc=isocolors(X,Y,Z,permute(mX,[2,1,3]),matsurf);
    nnz=nc==0;
    nc(nc~=0)=nc(nc~=0)-min(nc(nc~=0));
    nc(nc~=0)=(nc(nc~=0)/max(nc(nc~=0)))*64;
    nc(nnz)=0;
    nc(:)=ea_robustmean(round(nc)); % this will render them homogeneous.

    try
        jetlist=parula;
    catch
        jetlist=jet;
    end

    if exist('usecolormap','var')
        if ~isempty(usecolormap)
            jetlist=eval(usecolormap);
        end
    end

    jetlist=[0,0,0;jetlist];
    rgbnc=jetlist(round(nc)+1,:);
    set(matsurf,'FaceVertexCData',rgbnc);

    set(matsurf,'DiffuseStrength',0.9)
    set(matsurf,'SpecularStrength',0.1)
    set(matsurf,'FaceAlpha',0.2);
end

if showlabels
    for lab=1:length(atlasindices)
        [xx,yy,zz]=ind2sub(size(pV.img),find(pV.img==atlasindices(lab)));
        XYZ=[xx,yy,zz];
        centr_vx=[mean(XYZ,1),1]';
        centr_mm=pV.mat*centr_vx;

        labels(lab)=text(centr_mm(1),centr_mm(2),centr_mm(3),ea_underscore2space(atlaslegend{atlasindices(lab)}));
    end
end
