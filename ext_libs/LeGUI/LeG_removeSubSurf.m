function S = LeG_removeSubSurf(S)
%method to remove small surfaces enclosed within a main surrounding
%surface
%
%20210206
%TD

sv = S.vertices;
sf = S.faces;

minpts = size(sv,2)+1;

pd = pdist2(sv,sv,'euclidean','smallest',minpts+1);
epsln = ceil(prctile(pd(end,:),99));

idx = dbscan(sv,epsln,minpts); %epsilon is calculated by finding 99th percentile of distances for the minpts closest points
uidx = unique(idx);
len = nan(length(uidx),1);
for m=1:length(uidx)
    len(m) = sum(uidx(m)==idx);
end
[~,midx] = max(len);
ridx = uidx(setdiff(1:length(uidx),midx)); %indices to remove
ridx = find(ismember(idx,ridx));
kidx = setdiff(1:size(sv,1),ridx); %indices to keep
kidx_adj = 1:size(sv,1)-length(ridx); %indices to keep adjusted for shift

sv(ridx,:) = [];
sf(any(ismember(sf,ridx),2),:) = [];

[~,b] = ismember(sf(:),kidx);
sf = reshape(kidx_adj(b),[],3);

S.vertices = sv;
S.faces = sf;

% figure;
% pH = patch(S);
% set(pH,'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1,'AlphaDataMapping','none','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.3,'DiffuseStrength',0.8,'SpecularStrength',0.1,'SpecularExponent',10,'SpecularColorReflectance',1)
% 
% view([115,15])
% axis('vis3d','equal','tight')
% camlight;
% set(gca,'xgrid','on','ygrid','on','zgrid','on')
% set(gca,'xticklabelmode','auto','yticklabelmode','auto','zticklabelmode','auto')
% xlabel(gca,'x'); ylabel(gca,'y'); zlabel(gca,'z');
