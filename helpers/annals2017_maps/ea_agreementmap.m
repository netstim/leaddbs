function map=ea_agreementmap(niis,outputfilename)

for maps=1:length(niis)-1
    if ~exist('map','var')
        first=ea_load_nii(niis{maps});
    else
        first=map;
    end
    second=ea_load_nii(niis{maps+1});
    
    map=first;
    map.img((first.img.*second.img)<0)=nan; % set nonagreeing voxels to nan.
    map.img(map.img>0)=first.img(map.img>0).*second.img(map.img>0); % multiply positives
    map.img(map.img<0)=-first.img(map.img<0).*second.img(map.img<0); % multiply negatives
end
if exist('outputfilename','var')
    map.fname=outputfilename;
    ea_write_nii(map);
end

