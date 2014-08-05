function slicebw=ea_largestcomponent(slicebw)

try
stats=bwconncomp(slicebw);
catch
    disp('No image processing toolbox available. Using ls_bwconncomp by Stanis?aw Adaszewski.');
    try
stats=ls_bwconncomp(slicebw);
    


if stats.NumObjects>1
    maxlen=0;
    for i=1:stats.NumObjects
    if length(stats.PixelIdxList{i})>maxlen
        maxlen=length(stats.PixelIdxList{i});
        biggestobj=i;
    end
    end
    
    slicebw(:)=0;
    slicebw(stats.PixelIdxList{biggestobj})=1;
end
    


