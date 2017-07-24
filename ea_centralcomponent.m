function slicebw=ea_centralcomponent(slicebw,mask,options)

%% determine mask center:
[yy,xx]=find(mask);
maskcenter=[mean(xx),mean(yy)];

stats=ea_conncomp(slicebw);

if stats.NumObjects>1
    maxdist=10000;
    for i=1:stats.NumObjects
    if length(stats.PixelIdxList{i})>10
        
        sliceobj=slicebw;    
        sliceobj(:)=0;
        sliceobj(stats.PixelIdxList{i})=1;
        cen=ea_centroid(sliceobj);
        dist=ea_pdist([cen.Centroid;maskcenter]);
        ea_showdis(['Obj number ',num2str(i),', distance: ',num2str(dist),'.'],options.verbose);
        if dist<maxdist
        maxdist=dist;
        bestobj=i;
        ea_showdis(['Using object ',num2str(bestobj),'.'],options.verbose);
        
        end
    end
    end
    
    
    
    if exist('bestobj','var')
    
    slicebw(:)=0;
    slicebw(stats.PixelIdxList{bestobj})=1;
    else
        slicebw=ea_largestcomponent(slicebw);
    end
end
    

