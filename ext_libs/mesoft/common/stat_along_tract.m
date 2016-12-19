function bundle = stat_along_tract(ftr,mrdata)

if isempty(ftr),
    bundle.name = '';   
    bundle.fibcnt = 0;
    bundle.mean_length = nan;
    bundle.sdev_length  = nan;
    bundle.max_length = nan;
    bundle.min_length = nan;
    bundle.prctile05_length = nan;
    bundle.prctile95_length = nan;
    bundle.contrast = [];
    return
end;


for k = 1:length(ftr.fiber),
    fibs = ftr.curveSegCell(ftr.fiber{k}.curveID);
    fibs = cellfun(@(x) x*ftr.hMatrix(1:3,1:3)' + repmat(ftr.hMatrix(1:3,4)',[size(x,1) 1]),fibs,'uniformoutput',false); 
    fiblen = cellfun(@(x) sum(sqrt(sum((x(2:end,:)-x(1:end-1,:)).^2,2))),fibs);
    
    if isfield(ftr.fiber{k}.user,'overlap'),
        weight = overlap_to_weight(ftr.fiber{k}.user);
    else
        weight = ones(length(ftr.fiber{k}.curveID),1);
    end;
    
    bundle(k).name = ftr.fiber{k}.name;    
    bundle(k).fibcnt = sum(weight);
    bundle(k).mean_length = fiblen(:)'*weight(:) / sum(weight);
    bundle(k).sdev_length  = sqrt(((fiblen(:)- bundle(k).mean_length).^2)'*weight(:) / sum(weight));
    bundle(k).max_length = max(fiblen);
    bundle(k).min_length = min(fiblen);    
    bundle(k).prctile05_length = prctile(fiblen,5);
    bundle(k).prctile95_length = prctile(fiblen,95);
    bundle(k).contrast = [];
end;
    


if not(iscell(mrdata)),
    mrdata = {mrdata};
end;

for j = 1:length(mrdata),

    mr = mrdata{j};

    q = diag([-1 -1 1 1]);
    M = inv(mr.edges)*ftr.hMatrix;

    for k = 1:length(ftr.fiber),
        fibs = ftr.curveSegCell(ftr.fiber{k}.curveID);
        lens = cellfun(@(x) size(x,1),fibs);

    
        if isfield(ftr.fiber{k}.user,'overlap'),
            weight = overlap_to_weight(ftr.fiber{k}.user);
        else
            weight = ones(length(ftr.fiber{k}.curveID),1);
        end;
        weight = weight/sum(weight(:)) *length(ftr.fiber{k}.curveID);
            
        
        vertex = cat(1,fibs{:})-1;
        vertex(:,4) = 1; vertex = vertex*M'; vertex = vertex(:,1:3) +1;
        sz = size(mr.dataAy);
        vertex(vertex(:,1)<1,1) = 1;
        vertex(vertex(:,2)<1,2) = 1;
        vertex(vertex(:,3)<1,3) = 1;
        vertex(vertex(:,1)>sz(1),1) = sz(1);
        vertex(vertex(:,2)>sz(2),2) = sz(2);
        vertex(vertex(:,3)>sz(3),3) = sz(3);   
        vals = interp3(mr.dataAy,vertex(:,2),vertex(:,1),vertex(:,3));


        cumlens = [ 0 ;cumsum(lens)];

        N = 150;
        fi = zeros(N,length(fibs));
        for r = 1:length(fibs),
            fv{r} = vals(cumlens(r)+1:cumlens(r+1))'*weight(r);
            fi(:,r) = interp1(fv{r},1+(0:(N-1))/(N-1)*(length(fv{r})-1));
        end;       
        bundle(k).contrast(j).mean = mean(fi(:));
        bundle(k).contrast(j).sdev = std(fi(:));
        bundle(k).contrast(j).min = min(fi(:));
        bundle(k).contrast(j).max = max(fi(:));
        bundle(k).contrast(j).prctile05 = prctile(fi(:),5);
        bundle(k).contrast(j).prctile50 = prctile(fi(:),50);
        bundle(k).contrast(j).prctile95 = prctile(fi(:),95);


    end;
end;

return;



function w = overlap_to_weight(u)

w = prod(atan(u.overlap*10)/pi*2,2);


