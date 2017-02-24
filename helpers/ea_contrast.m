function slice=ea_contrast(slice,contrast,offset)
if ~exist('contrast','var')
    contrast=1;
end
if ~exist('offset','var')
    offset=0;
end

slice(slice(:)~=0)=contrast*zscore(slice(slice(:)~=0));
slice=slice+offset;
slice(slice>3)=3; % cut at 3 std devs
slice(slice<-3)=-3; % cut at -3 std devs

%slice=ea_sigmoid(slice);
slice=ea_minmax(slice);



function slice=ea_minmax(slice)
mn=min(slice(:));

slice=slice-mn;
mx=max(slice(:));
slice=slice/mx;

function g = ea_sigmoid(z)
g = 1.0 ./ (1.0 + exp(-z));
