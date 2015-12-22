function mstruc = createCM_GT(ftr,rois,Nsz,centroids)


display('building connect mat');
lens = single(cellfun(@len,ftr.curveSegCell));

fibs = ftr.curveSegCell;

terms = cellfun(@(x) [x(1,:); x(end,:)],fibs,'UniformOutput',false);
terms = cat(1,terms{:})';

rois(isnan(rois(:))) = 0;

if not(exist('centroids')),
    display('computing centroids');
    idx = setdiff(unique(rois(:)),0);
    [X Y Z] = ndgrid(1:size(rois,1),1:size(rois,2),1:size(rois,3));
    for k = 1:length(idx),
        indi = rois(:)==idx(k);
        indi = indi/sum(indi(:));
        centroids(idx(k),1) = sum(X(:).*indi(:));
        centroids(idx(k),2) = sum(Y(:).*indi(:));
        centroids(idx(k),3) = sum(Z(:).*indi(:));
    end;
end;

[cc lens] = CreateConnectivityMatrixROI(single(rois),double(Nsz),single(terms-1),single(lens));

fprintf('\n');
lens(lens>0) = lens(lens>0) ./ cc(lens>0);

mstruc.cc = cc;      %% connec matrix
mstruc.lens = lens;  %% lenght matrix
mstruc.centroids = centroids;

function l = len(x)

l = sum(sqrt(sum((x(2:end,:)-x(1:end-1,:)).^2,2)));













