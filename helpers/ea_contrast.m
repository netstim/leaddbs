function slice=ea_contrast(slice,contrast,offset)
% Enhance slice contrast

if ~exist('contrast','var')
    contrast=1;
end

if ~exist('offset','var')
    offset=0;
end

if ~isa(slice,'double')
    slice=double(slice);
end

ispositive=ea_nanmin(slice(:))>=0;
if sum(slice(:)~=0)
    slice(slice(:)~=0)=contrast*ea_nanzscore_sampled(slice(slice(:)~=0),15000);
end

if ispositive % only positive values
    slice(slice(:)==0)=ea_nanmin(slice(:));
end

slice=slice+offset;
slice(slice>3)=3; % cut at 3 std devs if above
slice(slice<-3)=-3; % cut at -3 std devs if above
slice=ea_rescale(slice);
slice=slice+offset;
slice=slice.*contrast;
