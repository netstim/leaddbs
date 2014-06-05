function [fitline,trajvector,diams]=reconstruct_trajectory(priorfitline,tra_nii,patientname,side,refine,options)


%Vtra=spm_vol([options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii']);
slice=zeros(size(tra_nii.img,1),size(tra_nii.img,2));
masknii=tra_nii;
if side==1
    if options.verbose>1; progressfig=figure('name','Finding left electrode'); end
else
    if options.verbose>1; progressfig=figure('name','Finding right electrode'); end
end

%maximize(progressfig)
if options.verbose>1; colormap(gray); end
endcount=0;
nanflag=0;

if ~refine % if this is not a refine-run but an initial run, mask of first slice has to be defined heuristically.
    % define initial mask
    mask=zeros(size(slice,1),size(slice,2));
    
    mask(200:350,70:220)=1;
    if side==2
        mask=fliplr(mask);
    end
    
    % initialize slice. mean average for entrypoint over the first 4 slices.
    slice=zeros(size(mask,1),size(mask,2),4);
    slicebw=zeros(size(mask,1),size(mask,2),4);
    for i=1:4
        [slice(:,:,i),slicebw(:,:,i)]=prepare_slice(tra_nii,mask,1,size(tra_nii.img,3)-(i-1),options);
    end
    slice=mean(slice,3);
    slicebw=logical(mean(slicebw,3));
    slicebw=centralcomponent(slicebw,mask,options);
    
    
    
    %keyboard % here to analyse initial slice.
    
    try
        stats=regionprops(slicebw,'Centroid');
    catch
        keyboard
    end
    try
        isempty(stats.Centroid); % this is only to check if stats.Centroid is empty.
        centerline(1,:)=[stats.Centroid,size(tra_nii.img,3)];
    catch
        disp('Threshold too high?');
        %pause;
    end
    
    
end





for sliceno=2:size(tra_nii.img,3) % sliceno is the counter (how many slices have been processed).
    
    if refine
        centerline(1,:)=priorfitline(1,:); % define initial point and mask for this run.
        mask=zeros(size(slice,1),size(slice,2));
        try
        estpoint=priorfitline(sliceno,:); % overwrite estpoint defined at the end of the loop if priorfitline is defined.
        catch
            break
        end
        mask(round(estpoint(2)-options.maskwindow):round(estpoint(2)+options.maskwindow),round(estpoint(1)-options.maskwindow):round(estpoint(1)+options.maskwindow))=1;
    end
    
    
    imgsliceno=size(tra_nii.img,3)-(sliceno-1); % imgsliceno is the slice number in the image.
    
    if imgsliceno<20
        showdis('Lower than slice 20. Stopping.',options.verbose);
        break
    end
    
    %% loop over each slice to find the electrode trajectory.
    %% part 1: finding the electrode on the current slice.
    %-------------------------------------------------------------------------------------------------%
    
    
    
    
    
    
    
    
    if options.verbose>1; subplot(3,3,1); end
    
    
    % the following function will return the slice and a bw copy of the
    % slice.
    [slice,slicebw,maskslice,maskslicebw]=prepare_slice(tra_nii,mask,sliceno,imgsliceno,options);
    
    if isempty(find(slicebw, 1)) % -> slice is not empty
    end
    
    
    
    % slice is always the raw current slice.
    
    if options.verbose>1; colormap(gray); end
    if options.verbose>1; imagesc(slice); end
    if options.verbose>1; colormap(gray); end
    if options.verbose>1; axis square; end
    
    
    
    
    
    %% part 2: check wheter the new point is plausible.
    %-------------------------------------------------------------------------------------------------%
    
    
    % the following function will use midpoint from the
    % last iteration to determine the distance to the new
    % one and will output the new midpoint.
    
    % check if estpoint has been defined
    
    
    
    
    
    
    if exist('estpoint','var')
        
        % this function will return one midpoint from the slice. If there are
        % more objects, it will return the midpoint of the one closest to the
        % estimated one.
        [numidpoint,diams(imgsliceno),greymaskslicebw]=findonemidpoint(slicebw,estpoint(1:2),mask,options);
        if isnan(numidpoint)
            showdis(['Midpoint is nan. Stopping.'],options.verbose);
            
            break
        end
        
        if pdist([estpoint;[numidpoint,imgsliceno]])<10-maxthree(refine)
            centerline(sliceno,:)=[numidpoint,imgsliceno];
            showdis(['Empirical Midpoint seems to be ',num2str([numidpoint,imgsliceno]),'.'],options.verbose);
            showdis(['New Midpoint found. Distance is ',num2str(pdist([estpoint;[numidpoint,imgsliceno]])),'.'],options.verbose);
        else
            
            endcount=endcount+1;
            if endcount==options.endtolerance
                showdis(['Too many interpolations. Stopping.'],options.verbose);
                
                break
            end
            centerline(sliceno,:)=estpoint;
            
            showdis(['No new Midpoint found. Distance is ',num2str(pdist([estpoint;[numidpoint,imgsliceno]])),'. Interpolating.'],options.verbose);
            
            if options.slow;            pause(1); end
        end
        
        showdis(['Diameter of this point is ',num2str(diams(imgsliceno)),'.'],options.verbose);
    else
        showdis('Estimated point not yet defined. Using second empirical point.',options.verbose);
        numidpoint=findonemidpoint(slicebw,centerline(1,1:2),mask,options);
        centerline(sliceno,:)=[numidpoint,imgsliceno];
    end
    
    
    
    
    
    
    
    
    
    
    
    %[centerline(sliceno,:),slice,nanflag,endslice]=determine_midpoint(tra_nii,centerline(sliceno-1,:),imgsliceno,sliceno,slice,slicebw,stats,centerline,options);
    
    
    if options.slow; pause(1); end
    
    
    
    
    
    
    %% part 3: update parameters for next run...
    %-------------------------------------------------------------------------------------------------%
    
    
    
    
    % this function estimates a fitted line and the following point based on the last points.
    
    [fitline,trajvector,estpoint]=fit_line(centerline);
    showdis(['Next point was estimated to be ',num2str(estpoint),'.'],options.verbose);
    % update mask
    mask=zeros(size(slice,1),size(slice,2));
    
    
    if round(estpoint(2)-options.maskwindow)<0 || round(estpoint(2)+options.maskwindow)>500 || round(estpoint(1)-options.maskwindow)<0 || round(estpoint(1)+options.maskwindow)>500
        display('Mask out of bounds. Must have lost trajectory. Try iterating with a different maskwindow again.');
        return
        %pause
    end
    
    
    mask(round(estpoint(2)-options.maskwindow):round(estpoint(2)+options.maskwindow),round(estpoint(1)-options.maskwindow):round(estpoint(1)+options.maskwindow))=1;
    
    
    
    %% part 4: visualization...
    %-------------------------------------------------------------------%
    
    if options.verbose>1; subplot(3,3,2); end
    if options.verbose>1; imagesc(mask); end
    if options.verbose>1; axis square; end
    
    if options.verbose>1; subplot(3,3,3); end
    if options.verbose>1; imagesc(slicebw); end
    if options.verbose>1; axis square; end
    
    
    if options.verbose>1; subplot(3,3,4); end
    
    if options.verbose>1; imagesc(maskslice); end
    if options.verbose>1; axis square; end
    
    if options.verbose>1; subplot(3,3,5); end
    if options.verbose>1; imagesc(maskslicebw); end
    if options.verbose>1; axis square; end
    
    
    if exist('greymaskslicebw','var')
        if ~isnan(greymaskslicebw)
            if options.verbose>1; subplot(3,3,6); end
            if options.verbose>1; imagesc(greymaskslicebw); end
            if options.verbose>1; axis square; end
        end
    else
        if options.verbose>1; subplot(3,3,6); end
        if options.verbose>1; imagesc(zeros(10)); end
        if options.verbose>1; axis square; end
    end
    
    
    
    
    
    
    
    
    if options.verbose>1; subplot(3,3,7:9); end
    
    
    if options.verbose>1; plot3(centerline(:,1),centerline(:,2),centerline(:,3)); hold on; plot3(centerline(1,1),centerline(1,2),centerline(1,3),'*g'); hold off; end
    
    
    if options.verbose>1; hold on; plot3(fitline(:,1),fitline(:,2),fitline(:,3),'r'); hold off; end
    
    if options.verbose>1; drawnow; end
    
    
    if options.slow; pause(0.1); end
    
    
    
end

if options.verbose>2; close(progressfig); end



function output=maxthree(refine) % simply returns maximally four
output=refine;
if output>3; output=3; end