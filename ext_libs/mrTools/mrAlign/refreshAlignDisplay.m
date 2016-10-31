function img = refreshAlignDisplay(handles)

%   $Id$

global ALIGN

% Just return if volume hasn't been loaded yet.
if isempty(ALIGN.volume)
    return
end
    
% Make mrAlign the current figure and raise it
figure(handles.figure1)

% Compose xform
xform = ALIGN.guiXform * ALIGN.xform;

% Extract base and overlay images and convert to RGB
sliceOrientation = ALIGN.sliceOrientation;
sliceNum = ALIGN.coords(sliceOrientation);
alpha = get(handles.transparencySlider,'Value');
if ~isempty(ALIGN.volume)
    switch sliceOrientation
        case 1
            base = squeeze(ALIGN.volume(sliceNum,:,:));
        case 2
            base = squeeze(ALIGN.volume(:,sliceNum,:));
        case 3
            base = squeeze(ALIGN.volume(:,:,sliceNum));
    end
    baseRGB = rescale2rgb(base,ALIGN.baseCmap,ALIGN.volumeClip);
end
if ~isempty(ALIGN.inplanes)
    overlay = transformInplanes(ALIGN.inplanes,ALIGN.volSize,xform,...
        sliceOrientation,sliceNum);
    overlayRGB = rescale2rgb(overlay,ALIGN.overlayCmap,ALIGN.inplanesClip);
end

% Figure out what to display depending on status of overlay button and
% whether or not base and overlay images exist.
if get(handles.overlayButton,'Value')
	if exist('baseRGB','var') & exist('overlayRGB','var')
		alphaMap = repmat(alpha * ~isnan(overlay),[1 1 3]);
		img = (1-alphaMap).*baseRGB + alphaMap.*overlayRGB;
	elseif exist('baseRGB','var')
		img = baseRGB;
	elseif exist('overlayRGB','var')
		img = overlayRGB;
	end
elseif exist('baseRGB','var')
	img = baseRGB;
end

% Transpose and flip
if exist('img','var')
    if get(handles.transposeButton,'Value')
        imgSize = size(img);
        imgTranspose = zeros(imgSize(2), imgSize(1), imgSize(3));
        for c = 1:3
            imgTranspose(:,:,c) = img(:,:,c)';
        end
        img = imgTranspose;
    end
    if get(handles.flipButton,'Value')
        img = flipdim(img,1);
	end
else
	img = [];
end

% Display it
if ~isempty(img)
    image(img);
    axis off
    axis image
end

% *** Something like this should work (but doesn't) using matlab's built-in
% transparency functionality. So I implemented the alphamap from scratch
% above.
%     baseImage = subimage(base,ALIGN.baseCmap);
%     hold on
%     overlayImage = subimage(overlay,ALIGN.overlayCmap);
%     set(overlayImage,'AlphaData',alpha);
%     hold off
%     axis off



function rgb = rescale2rgb(image,cmap,clip)
%function rgb = rescale2rgb(image,cmap,[clipMin,clipMax])

% Clip
clipMin = clip(1);
clipMax = clip(2);
image(image < clipMin) = clipMin;
image(image > clipMax) = clipMax;

% Scale
indices = round(255 * (image-clipMin)/(clipMax-clipMin)) + 1;
indices = max(1,min(indices,size(cmap,1)));

% Extract r,g,b components
r = zeros(size(image));
g = zeros(size(image));
b = zeros(size(image));
r(:) = cmap(indices,1);
g(:) = cmap(indices,2);
b(:) = cmap(indices,3);

% Stuff them into rgb
dims = [size(image),3];
rgb = cat(length(dims),r,g,b);


function inplanesXform = transformInplanes(inplanes,volSize,xform,sliceOrientation,slice)

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;
xform = inv(shiftXform) * inv(xform) * shiftXform;
%let's round the xform matrix for display. This avoids missing voxels
%when the transforms of source and destination should be identical but their 
%division results in non-zero rotations
xform = 1e-7*round(1e7*xform);

% Generate coordinates with meshgrid
switch sliceOrientation
    case 1
        x = slice * ones(volSize(2),volSize(3));
        [z,y] = meshgrid(1:volSize(3),1:volSize(2));
    case 2
        y = slice * ones(volSize(1),volSize(3));
        [z,x] = meshgrid(1:volSize(3),1:volSize(1));
    case 3
        z = slice * ones(volSize(1),volSize(2));
        [y,x] = meshgrid(1:volSize(2),1:volSize(1));
end
dims = size(x);
numPixels = prod(dims);
xvec = reshape(x,1,numPixels);
yvec = reshape(y,1,numPixels);
zvec = reshape(z,1,numPixels);
coords = [xvec; yvec; zvec; ones(1,numPixels)];

% Transform coordinates
coordsXform = xform * coords;
xi = reshape(coordsXform(1,:),dims);
yi = reshape(coordsXform(2,:),dims);
zi = reshape(coordsXform(3,:),dims);

% Interpolate
% Note: interp3 treats x and y in right-handed coordinate system, not in
% matrix index order so we need to swap them here. See example code
% interpVolume.m.
inplanesXform = interp3(inplanes,yi,xi,zi,'linear',nan);
