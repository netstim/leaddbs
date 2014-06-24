function [slice,slicebw,maskslice,maskslicebw]=ea_prepare_slice(nii,mask,sliceno,imgsliceno,options)
% this function binarizes each axial slice (resulting in plane slicebw) by
% thresholding the variable slice by x times the standard-deviation of its
% intensity values. x here refers to the sum of a factor defined in the
% options of the program (options.tra_stdfactor) and a value adjusting it
% if it doesn't seem to be properly set for the image.

slice=double(nii.img(:,:,imgsliceno))'; % extract the correct slice.

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


threshold=mean(maskslice(:))-(options.tra_stdfactor+add_to_tra_stdfactor)*std(maskslice(:)); % =80.


%% make binary thresholded copies of slice.
slicebw=zeros(length(slice));
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


threshold=nanmean(maskslice(:))-(options.tra_stdfactor+addfactor)*nanstd(maskslice(:)); % =80.


%% make binary thresholded copies of slice.
slicebw=zeros(length(slice));
slicebw(slice<threshold)=1;
slicebw(~mask)=0;



if isempty(find(slicebw(:), 1))
   addfactor=addfactor-cnt*0.1;
   ea_showdis(['Lowering initial factor to ',num2str(addfactor),'.'],options.verbose);
    
else
    found=1;
end
    

end

