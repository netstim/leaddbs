function map=ea_agreementmap(nii1,nii2,outputfilename);

first=ea_load_nii(nii1);
second=ea_load_nii(nii2);

map=first;
map.img((first.img.*second.img)<0)=nan; % set nonagreeing voxels to nan.
map.img(map.img>0)=first.img(map.img>0).*second.img(map.img>0); % multiply positives
map.img(map.img<0)=-first.img(map.img<0).*second.img(map.img<0); % multiply negatives

map.fname=outputfilename;
ea_write_nii(map);

