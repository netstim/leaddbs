function [trajectory,trajvector]=ea_reconstruct_trajectory(priortrajectory,tra_nii,side,refine,options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if strcmp(options.subj.postopModality, 'CT') % CT support
    tra_nii.img=tra_nii.img*-1;
else
    if strcmp(options.entrypoint,'Auto')
        warning('Automatic entry point detection not implemented for MRI. Setting to Manual.')
        options.entrypoint = 'Manual';
    end
end

%Vtra=spm_vol([options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii']);
slice=zeros(size(tra_nii.img,2),size(tra_nii.img,1));
masknii=tra_nii;

if options.verbose>1
    colormap(gray);
end
endcount=0;
nanflag=0;

% determine startslice at ~ z=8.7mm
[startslice,endslice,masksz] = ea_getstartslice(options);
mmpt = [0;0;startslice;1];
mmvx = tra_nii.mat\mmpt;
startslice = round(mmvx(3));
clear mmpt mmvx

if side > 2 % always go for manual entrypoint in case of >second electrode.
    options.entrypoint='Manual';
end

flipside = 1+(tra_nii.mat(1)<0);

if ~refine % if this is not a refine-run but an initial run, mask of first slice has to be defined heuristically.
    % define initial mask
    mask = zeros(size(slice,1),size(slice,2));
    switch options.entrypoint
        case 'Manual'
            colormask = zeros(size(slice,1),size(slice,2),3);
            colormask(:,:,1) = 1;
            mask(masksz(1):masksz(2),masksz(3):masksz(4)) = 1;

            if side == flipside
                mask = flip(mask,2);
            end

            slice = double(tra_nii.img(:,:,startslice))'; % extract the correct slice.
            %slice=fliplr(slice);
            slice(slice==0) = nan;
            mn = figure('color','w','ToolBar','none','NumberTitle','off','Menubar','none','name','Please specify manual starting point.');
            ax = axes;
            WinOnTop(mn, true); % bring window to top
            set(0, 'CurrentFigure', mn);    % set current figure explicitly
            set(mn, 'CurrentAxes', ax);    % set current axes explicitly
            imagesc(slice);
            colormap(gray);
            hold on
            cof = imshow(colormask);
            set(cof, 'AlphaData', mask*0.3);
            hold off
            [X, Y] = ginput(1);
            close(mn);
            % reset mask from mouse input
            mask = zeros(size(slice,1),size(slice,2));
            mask(round(Y-10:Y+10),round(X-10:X+10)) = 1;
        case 'Auto'
            % TP: Similar to manual, but try to detect entry-point using artefact
            crop_n = 100;
            slice = double(tra_nii.img(:,:,startslice))'; % extract the correct slice.
            [h,w] = size(slice);
            crop = slice(crop_n:h-crop_n, crop_n:w-crop_n);
            midpt = floor(size(crop,2)/2);
            idx = 1:midpt;
            if side == flipside
                idx = midpt:size(crop,2);
            end
            b = crop(:,idx);
            s = b;
            s(s>(min(b(:))*0.2)) = 0;
            height = abs(min(s(:))) * 0.05;
            [~, xid] = findpeaks(-sgolayfilt(min(s,[],1), 1, 21), 'MinPeakHeight', height);
            [~, yid] = findpeaks(-sgolayfilt(min(s,[],2), 1, 21), 'MinPeakHeight', height);
            if (idx(1) == 1)
                [~, id] = min(abs(xid-size(crop,2)));
            else
                [~, id] = min(abs(xid));
            end
            X = xid(id);
            [~, id] = min(abs(yid-size(crop,1)));
            Y = yid(id) + crop_n;
            X = X + idx(1)-1 + crop_n;

            % further refine
            crop_n = 20;
            crop = slice(Y-crop_n:Y+crop_n, X-crop_n:X+crop_n);
            [~, idx] = min(crop(:));
            [Y2, X2] = ind2sub(size(crop), idx);
            X = round(X + X2 - size(crop,2)/2);
            Y = round(Y + Y2 - size(crop,1)/2);

            mask = zeros(h,w);
            mask(round(Y-10:Y+10),round(X-10:X+10)) = 1;
        otherwise
            mask(masksz(1):masksz(2),masksz(3):masksz(4)) = 1;
            if side == flipside
                mask = fliplr(mask);
            end
    end

    % initialize slice. mean average for entrypoint over the first 4 slices.
    slice = zeros(size(mask,1),size(mask,2),4);
    slicebw = zeros(size(mask,1),size(mask,2),4);

    for i = 10:14
        try
        [slice(:,:,i),slicebw(:,:,i)] = ea_prepare_slice(tra_nii,mask,1,startslice-(i-1),options);
        catch
            % keyboard % TP: Requesting input here pauses tasks when batch
            % processing. I made a work around,but not sure if it will lead
            % to further errors.
            if (length(options.uipatdirs) > 1)
                warning([options.patientname, ': Could not prepare slice when reconstructing trajectory.']);
            else
                keyboard
            end
        end
    end
    slice = mean(slice,3);
    slicebw = logical(mean(slicebw,3));
    slicebw = ea_centralcomponent(slicebw,mask,options);

    %keyboard % here to analyse initial slice.

    stats = ea_centroid(slicebw);
    try
        isempty(stats.Centroid); % this is only to check if stats.Centroid is empty.
        centerline(1,:) = [stats.Centroid,startslice];
    catch
        disp('Threshold too high?');
    end
end

Vmat = nii2Vmat(tra_nii);
zfifteen = Vmat\[0;0;endslice;1];

if side == 1
    if options.verbose>1
        progressfig = figure('name','Finding right electrode','NumberTitle','off','Menubar','none','ToolBar','none');
        set(gcf,'color','w');
        axis off;
    end
elseif side == 2
    if options.verbose>1
        progressfig = figure('name','Finding left electrode','NumberTitle','off','Menubar','none','ToolBar','none');
        set(gcf,'color','w');
        axis off;
    end
else
    if options.verbose > 1
        progressfig = figure('name',['Finding electrode number ',num2str(side)],'NumberTitle','off','Menubar','none','ToolBar','none');
        set(gcf,'color','w');
        axis off;
    end
end

set(progressfig,'KeyPressFcn',@ea_keystr);

%% starting slice 2:end
for sliceno = 2:startslice % sliceno is the counter (how many slices have been processed).
    % uncomment the following two lines to write out slice views.
    %imwrite(((reshape(slice(logical(mask)),sqrt(numel(find(mask))),sqrt(numel(find(mask)))))-min(slice(:)))/(max(slice(:)-min(slice(:)))),['slice_',num2str(sliceno),'.png']);
    %imwrite(reshape(slicebw(logical(mask)),sqrt(numel(find(mask))),sqrt(numel(find(mask)))),['slicebw_',num2str(sliceno),'.png']);
    if refine
        centerline(1,:) = priortrajectory(1,:); % define initial point and mask for this run.
        mask=zeros(size(slice,1),size(slice,2));
        try
            estpoint = priortrajectory(sliceno,:); % overwrite estpoint defined at the end of the loop if priortrajectory is defined.
        catch
            break
        end
        mask(round(estpoint(2)-options.maskwindow):round(estpoint(2)+options.maskwindow),round(estpoint(1)-options.maskwindow):round(estpoint(1)+options.maskwindow))=1;
    end

    imgsliceno = startslice-(sliceno-1); % imgsliceno is the slice number in the image.
    if imgsliceno < zfifteen(3) && ~strcmp(options.entrypoint,'Cg25')
        ea_showdis('Lower than z=-15.5 mm. Stopping.',options.verbose);
        break
    end

    %% loop over each slice to find the electrode trajectory.
    %% part 1: finding the electrode on the current slice.
    %-------------------------------------------------------------------------------------------------%
    % the following function will return the slice and a bw copy of the
    % slice.
    [slice,slicebw,maskslice,maskslicebw] = ea_prepare_slice(tra_nii,mask,sliceno,imgsliceno,options);
    if isempty(find(slicebw, 1)) % -> slice is not empty
        % HANDLE EMPTY SLICE
    end

    % slice is always the raw current slice.
    if options.verbose>1
        ea_setfocus(progressfig);
        subplot(3,3,1);
        imagesc(slice);
        colormap(gray);
        axis off square;
    end

    %% part 2: check whether the new point is plausible.
    %-------------------------------------------------------------------------------------------------%
    % the following function will use midpoint from the
    % last iteration to determine the distance to the new
    % one and will output the new midpoint.

    % check if estpoint has been defined
    if exist('estpoint','var')
        % this function will return one midpoint from the slice. If there are
        % more objects, it will return the midpoint of the one closest to the
        % estimated one.
        [numidpoint,greymaskslicebw,options]=ea_findonemidpoint(slicebw,estpoint(1:2),mask,options);
        if isnan(numidpoint)
%            ea_showdis('Midpoint is nan. Stopping.', options.verbose);
%            break
            numidpoint=estpoint(1:2);
        end

        if ea_pdist([estpoint;[numidpoint,imgsliceno]])<15-maxthree(refine)
            centerline(sliceno,:)=[numidpoint,imgsliceno];
            %ea_showdis(['Empirical Midpoint seems to be ',num2str([numidpoint,imgsliceno]),'.'],options.verbose);
            %ea_showdis(['New Midpoint found. Distance is ',num2str(ea_pdist([estpoint;[numidpoint,imgsliceno]])),'.'],options.verbose);
        else
            endcount=endcount+1;
            if endcount==options.endtolerance
                ea_showdis('Too many interpolations. Stopping.',options.verbose);
                break
            end
            centerline(sliceno,:)=estpoint;
            %ea_showdis(['No new Midpoint found. Distance is ',num2str(ea_pdist([estpoint;[numidpoint,imgsliceno]])),'. Interpolating.'],options.verbose);
        end
    else
        ea_showdis('Estimated point not yet defined. Using second empirical point.',options.verbose);
        numidpoint=ea_findonemidpoint(slicebw,centerline(1,1:2),mask,options);

        centerline(sliceno,:)=[numidpoint,imgsliceno];
        if isnan(centerline)
            ea_error('Reconstruction failed. Please choose "manual" entrypoint.');
        end
    end

    endnow=getappdata(progressfig,'endnow');
    if ~isempty(endnow)
        if endnow
            ea_showdis('User pressed space, stopping.',options.verbose);
            break
        end
    end

    %% part 3: update parameters for next run...
    %-------------------------------------------------------------------------------------------------%
    % this function estimates a fitted line and the following point based on the last points.
    [trajectory,trajvector,estpoint]=ea_fit_line(centerline);

    % ea_showdis(['Next point was estimated to be ',num2str(estpoint),'.'],options.verbose);
    % update mask
    mask=zeros(size(slice,1),size(slice,2));
    if round(estpoint(2)-options.maskwindow) < 0 || ...
            round(estpoint(2)+options.maskwindow) > 500 || ...
            round(estpoint(2)+options.maskwindow) > size(slice,1) || ...
            round(estpoint(1)-options.maskwindow) < 0 || ...
            round(estpoint(1)+options.maskwindow) > 500 || ...
            round(estpoint(1)+options.maskwindow) > size(slice,2)
        close(progressfig)
        ea_error(sprintf(['Mask out of bounds! Must have lost trajectory...\n' ...
            'Please try a different ''Mask window size'' or ' ...
            'try manual mode by setting ''Entrypoint for Target'' to ''Manual''.\n'...
            'Patient: ', options.patientname]), ...
            'Electrode Reconstruction Error', ...
            dbstack);
        return
        %pause
    end
    % TP: Fixing bug. Sometimes when the trajectory is not found the index
    % to mask is less than 1 causing Lead to crash
    try
        mask(round(estpoint(2)-options.maskwindow : estpoint(2)+options.maskwindow), ...
            round(estpoint(1)-options.maskwindow : estpoint(1)+options.maskwindow))=1;
    catch ME
        if (strcmp(ME.identifier, 'MATLAB:badsubscript'))
            ea_error(sprintf(['Mask index out of bounds! Must have lost trajectory...\n' ...
            'Please try a different ''Mask window size'' or ' ...
            'try manual mode by setting ''Entrypoint for Target'' to ''Manual''.\n'...
            'Patient: ', options.patientname]), ...
            'Electrode Reconstruction Error', ...
            dbstack);
            return;
        end
    end
    %% part 4: visualization...
    %-------------------------------------------------------------------%
    if options.verbose>1
        ea_setfocus(progressfig);
        subplot(3,3,2);
        imagesc(mask);
        axis off square;

        subplot(3,3,3);
        imagesc(slicebw);
        axis off square;

        subplot(3,3,4);
        imagesc(maskslice);
        axis off square;

        subplot(3,3,5);
        imagesc(maskslicebw);
        axis off square;
    end

    if exist('greymaskslicebw','var')
        if ~isnan(greymaskslicebw)
            if options.verbose>1
                ea_setfocus(progressfig);
                subplot(3,3,6);
                imagesc(greymaskslicebw);
                axis off square;
            end
        end
    else
        if options.verbose>1
            ea_setfocus(progressfig);
            subplot(3,3,6);
            imagesc(zeros(10));
            axis off square;
        end
    end

    if options.verbose>1
        ea_setfocus(progressfig);
        subplot(3,3,7:9);
        plot3(centerline(:,1),centerline(:,2),centerline(:,3));
        hold on;
        plot3(centerline(1,1),centerline(1,2),centerline(1,3),'*g');
        plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3),'r');
        hold off;
        drawnow;
    end
end

if options.verbose>2
    close(progressfig);
end


function ea_setfocus(progressfig)
try
    set(0,'CurrentFigure',progressfig)
catch
    ea_error('Please do not close the progress figure during reconstruction.')
end


function ea_keystr(progressfig,event)
commnd=event.Character;
switch lower(commnd)
    case ' '
        setappdata(progressfig,'endnow',1);
end


function Vmat=nii2Vmat(nii)
Vmat=zeros(4);
Vmat(1,:)=nii.hdr.hist.srow_x-[0,0,0,nii.hdr.dime.pixdim(2)];
Vmat(2,:)=nii.hdr.hist.srow_y-[0,0,0,nii.hdr.dime.pixdim(3)];
Vmat(3,:)=nii.hdr.hist.srow_z-[0,0,0,nii.hdr.dime.pixdim(4)];
Vmat(4,:)=[0,0,0,1];


function output=maxthree(refine) % simply returns maximally four
output=refine;
if output>3
    output=3;
end


function [startslice,endslice,masksz]=ea_getstartslice(options) % get reconstruction default dimensions for current space
spacedef=ea_getspacedef;
standardspacedef=load([ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym',filesep,'spacedef.mat']);
if isfield(spacedef,'guidef')
    whichentry=ismember(options.entrypoint,spacedef.guidef.entrypoints);
    masksz=spacedef.guidef.masks(whichentry,:);
    startslice=spacedef.guidef.startslice;
    endslice=spacedef.guidef.endslice;
else % use MNI defaults
    startslice=8.7; % default height at where to start auto reconstruction
    endslice=-15.5;
    switch options.entrypoint
        case 'STN, GPi or ViM'
            masksz=standardspacedef.spacedef.guidef.masks(1,:);
        case 'Cg25'
            masksz=standardspacedef.spacedef.guidef.masks(2,:);
    end
end

if strcmp(options.entrypoint,'Manual') || strcmp(options.entrypoint,'Auto')
    try
        masksz=spacedef.guidef.masks(1,:);
    catch
        masksz=standardspacedef.spacedef.guidef.masks(1,:); % use STN default
    end
end

if ~exist('masksz','var') % e.g. in case manual set
    try
        masksz=spacedef.guidef.masks(1,:); % use first entry if defined
    catch
        masksz=standardspacedef.spacedef.guidef.masks(1,:); % use STN default
    end
end
