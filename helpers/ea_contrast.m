function slice=ea_contrast(slice,contrast,offset)
if ~exist('contrast','var')
    contrast=1;
end
if ~exist('offset','var')
    offset=0;
end
if ~strcmp(class(slice),'double')
    slice=double(slice);
end
%disp([num2str(contrast),',',num2str(offset)]);
ispositive=ea_nanmin(slice(:))>=0;
slice(slice(:)~=0)=contrast*ea_nanzscore_sampled(slice(slice(:)~=0),15000);

if ispositive; % only positive values
slice(slice(:)==0)=ea_nanmin(slice(:));
end
slice=slice+offset;
slice(slice>3)=3; % cut at 3 std devs if above
slice(slice<-3)=-3; % cut at -3 std devs if above
slice=ea_minmax(slice);
slice=slice+offset;
slice=slice.*contrast;




function slice=ea_minmax(slice)
mn=min(slice(:));

slice=slice-mn;
mx=max(slice(:));
slice=slice/mx;

function g = ea_sigmoid(z)
g = 1.0 ./ (1.0 + exp(-z));
