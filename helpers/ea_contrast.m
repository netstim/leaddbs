function slice=ea_contrast(slice,contrast,offset)


slice(:)=contrast*zscore(slice(:));
slice=slice+offset;

slice=ea_sigmoid(slice);




function slice=ea_minmax(slice)
mn=min(slice(:));

slice=slice-mn;
mx=max(slice(:));
slice=slice/mx;

function g = ea_sigmoid(z)
g = 1.0 ./ (1.0 + exp(-z));
