function [slice,slicebw,maskslice,maskslicebw]=ea_prepare_slice(nii,mask,sliceno,imgsliceno,options)
% this function binarizes each axial slice (resulting in plane slicebw) by
% thresholding the variable slice by x times the standard-deviation of its
% intensity values. x here refers to the sum of a factor defined in the
% options of the program (options.tra_stdfactor) and a value adjusting it
% if it doesn't seem to be properly set for the image.

slice=double(nii.img(:,:,imgsliceno))'; % extract the correct slice.
%slice=ea_smooth2(slice,5,5);
for lower=1:20
    if ~any(slice(:))
        slice=double(nii.img(:,:,imgsliceno-lower))';
    else
        break
    end
end

%slice=fliplr(slice);
slice(slice==0)=nan;
if sliceno<3 % first slices: lower tra_factor since great mask is used.
    add_to_tra_stdfactor=ea_determine_addfactor(1,slice,mask,options);
else
    add_to_tra_stdfactor=0;
end

maskslice=slice(logical(mask));

maskslice=reshape(maskslice,sqrt(length(maskslice)),sqrt(length(maskslice)));
maskslice=ea_smooth2(maskslice,3,3);

%figure, imagesc(maskslice)

threshold=mean(maskslice(:))-(options.tra_stdfactor+add_to_tra_stdfactor)*std(maskslice(:)); % =80.

%% make binary thresholded copies of slice.
slicebw=slice;
slicebw(:)=0;
slicebw(slice<threshold)=1;
slicebw(~mask)=0;

maskslicebw=zeros(length(maskslice));
maskslicebw(maskslice<threshold)=1;


function addfactor=ea_determine_addfactor(addfactor,slice,mask,options)

cnt=1;
found=0;

while ~found
    maskslice=slice(logical(mask));
    maskslice=reshape(maskslice,sqrt(length(maskslice)),sqrt(length(maskslice)));

    threshold=ea_nanmean(maskslice(:))-(options.tra_stdfactor+addfactor)*ea_nanstd(maskslice(:)); % =80.

    %% make binary thresholded copies of slice.
    slicebw=zeros(length(slice));
    slicebw(slice<threshold)=1;
    slicebw(~mask)=0;

    if isempty(find(slicebw(:), 1))
        addfactor=addfactor-cnt*0.1;
        ea_showdis(['Lowering initial factor to ',num2str(addfactor),'.'],options.verbose);
        if cnt>500
            ea_error('Trajectory could not be found. Please choose manual entry-point to select trajectory manually and make sure that images are properly normalized in MNI space.');
        end
    else
        found=1;
    end
end

function y = ea_nanmean(varargin)
if nargin==2
    x=varargin{1};
    dim=varargin{2};
elseif nargin==1
    x=varargin{1};
    dim=1;
end

N = sum(~isnan(x), dim);
y = ea_nansum(x, dim) ./ N;
