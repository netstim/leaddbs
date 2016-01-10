function [matsurf]=ea_showconnectivitypatch(resultfig,pV,mX,thresh)

tmX=mX>thresh;


bb=[0,0,0;size(mX)];

bb=map_coords_proxy(bb,pV);
gv=cell(3,1);
for dim=1:3
    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(mX,dim));
end
[X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
fv=isosurface(X,Y,Z,permute(tmX,[2,1,3]),0.5); % graph_metric
set(0,'CurrentFigure',resultfig)
matsurf=patch(fv,'facealpha',0.7,'EdgeColor','none','facelighting','phong','FaceColor','interp');

 zidx=mX==0;
 zidx=logical(zidx+isnan(mX));
% 
% %cgX(isnan(cgX))=0;
% try % if cgX isn't empty..
% cgX(cgX~=0)=cgX(cgX~=0)-ea_nanmin(cgX(cgX~=0));
% 
% cgX(cgX~=0)=(cgX(cgX~=0)/ea_nanmax(cgX(cgX~=0)))*64;
% end
% 
 mX(zidx)=0; % reset prior zero/nan values to zero

nc=isocolors(X,Y,Z,permute(mX,[2,1,3]),matsurf);
nnz=nc==0;
nc(nc~=0)=nc(nc~=0)-min(nc(nc~=0));
nc(nc~=0)=(nc(nc~=0)/max(nc(nc~=0)))*64;
nc(nnz)=0;


try
jetlist=parula;
catch
    jetlist=jet;
end
jetlist=[0,0,0;jetlist];
rgbnc=jetlist(round(nc)+1,:);
set(matsurf,'FaceVertexCData',rgbnc);

set(matsurf,'DiffuseStrength',0.9)
set(matsurf,'SpecularStrength',0.1)
set(matsurf,'FaceAlpha',0.2);














% %% old
% 
% bb=[0,0,0;size(mX)];
% 
% bb=map_coords_proxy(bb,pV);
% gv=cell(3,1);
% for dim=1:3
%     gv{dim}=linspace(bb(1,dim),bb(2,dim),size(mX,dim));
% end
% [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
% 
% 
% fv=isosurface(X,Y,Z,permute(mX,[2,1,3]),thresh); % connected regions
% 
% if ischar(options.prefs.lhullsimplify)
%     
%     % get to certain number of faces
%     simplify=10000/length(fv.faces);
%     fv=reducepatch(fv,simplify);
%     
% else
%     if options.prefs.lhullsimplify<1 && options.prefs.lhullsimplify>0
%         fv=reducepatch(fv,options.prefs.lhullsimplify);
%     end
% end
% 
% set(0,'CurrentFigure',resultfig)
% matsurf=patch(fv,'FaceColor','interp','facealpha',0.7,'EdgeColor','none','facelighting','phong');
% 
% cmX=mX;
% zidx=cmX==0;
% zidx=logical(zidx+isnan(cmX));
% cmX(isnan(cmX))=0;
% 
% try % if not possible, connectivity is empty..
% cmX=cmX-min(cmX(cmX~=0));
% end
% try % if not possible, connectivity is empty..
% cmX=(cmX/max(cmX(cmX~=0)))*64;
% end
% cmX(zidx)=0; % reset prior zero/nan values to zero
% isocolors(X,Y,Z,permute(cmX,[2,1,3]),matsurf)
% 
% set(matsurf,'DiffuseStrength',0.9)
% set(matsurf,'SpecularStrength',0.1)
% set(matsurf,'FaceAlpha',0.3);
% %threshold seedcon
% % seedcon=((seedcon-thresh)/(max(seedcon)-thresh))*255;
% % 
% % cX=pX;
% % for i=1:length(seedcon)
% %     cX(pX==i)=seedcon(i);
% % end
% %nc=isonormals(X,Y,Z,permute(mX,[2,1,3]),matsurf);
% %nc=isocolors(X,Y,Z,permute(cX,[2,1,3]),matsurf);
% %matsurf.FaceColor='interp';


function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';

