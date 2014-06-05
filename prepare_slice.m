function [slice,slicebw,maskslice,maskslicebw]=prepare_slice(nii,mask,sliceno,imgsliceno,options)
% this function binarizes each axial slice (resulting in plane slicebw) by
% thresholding the variable slice by x times the standard-deviation of its
% intensity values. x here refers to the sum of a factor defined in the
% options of the program (options.tra_stdfactor) and a value adjusting it
% if it doesn't seem to be properly set for the image.

slice=double(nii.img(:,:,imgsliceno))'; % extract the correct slice.
slice=fliplr(slice);
slice(slice==0)=nan;

if sliceno<3 % first slices: lower tra_factor since great mask is used.
add_to_tra_stdfactor=determine_addfactor(1,slice,mask,options);

else
add_to_tra_stdfactor=0;
end

%% adapt the threshold for the first 10 slices.
% if sliceno<10
%     % in first slices, check how many pixels were selected and if not, add a value to the tra_std_factor to select more pixels.
%     showdis('First slices, adapting threshold...',options.verbose);
% %     goodthresholdfound=0; % reset flag.
% %     while ~goodthresholdfound % flag will be set in adapt_threshold.
% %         % the following function corrects the threshold by
% %         % checking how many pixels have been found in this slice
% %         % (and wether this is enough).
% %         [add_to_tra_stdfactor,goodthresholdfound]=adapt_threshold(add_to_tra_stdfactor,slice,nii,mask,sliceno,options);
% %         
% %     end
%     
% end


    maskslice=slice(logical(mask));
    maskslice=reshape(maskslice,sqrt(length(maskslice)),sqrt(length(maskslice)));


threshold=mean(maskslice(:))-(options.tra_stdfactor+add_to_tra_stdfactor)*std(maskslice(:)); % =80.


%% make binary thresholded copies of slice.
slicebw=zeros(length(slice));
slicebw(slice<threshold)=1;


slicebw(~mask)=0;

maskslicebw=zeros(length(maskslice));
maskslicebw(maskslice<threshold)=1;

