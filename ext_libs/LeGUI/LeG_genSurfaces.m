function [BrainSurf,ProjSurf] = LeG_genSurfaces(ImageCell,ImgInfo) 
%Generates brain surface and a surface for projecting and computing a normal vector for moving an
%electrode. Two volumes are summed together (MRGray, MRWhite) and then thresholded. All volumes must
%be the same dimension and contain values from 0 to 1. ImageCell is a cell of the volumes and
%ImgInfo is the nifti header for one of the volumes.
%
%20200930

% XYZScale = [sqrt(sum(ImgInfo.mat(:,2).^2)),sqrt(sum(ImgInfo.mat(:,1).^2)),sqrt(sum(ImgInfo.mat(:,3).^2))]; %x/y reversed in spm world (this might be wrong?)
XYZScale = sqrt(sum(ImgInfo.mat(1:3,1:3).^2));

Img = zeros(size(ImageCell{1}));
for k=1:length(ImageCell)
    Img = Img + ImageCell{k};
end

Img = Img>0.95;

filtsizemm = 1.5;
filtsizevox = filtsizemm*ones(1,3).*1./XYZScale;
idx = mod(filtsizevox,2)<1;
filtsizevox = floor(filtsizevox);
filtsizevox(idx) = filtsizevox(idx)+1;
se = strel('cuboid',filtsizevox);

Img = imerode(Img,se);

CC = bwconncomp(Img);
S = regionprops(CC,'area');

A = [S.Area];
[~,idx] = max(A);
px = CC.PixelIdxList(setdiff(1:length(A),idx));

for k=1:length(px)
    Img(px{k}) = false;
end

Img = imdilate(Img,se);

filtsizemm = 3;
filtsizevox = filtsizemm*ones(1,3).*1./XYZScale;
idx = mod(filtsizevox,2)<1;
filtsizevox = floor(filtsizevox);
filtsizevox(idx) = filtsizevox(idx)+1;

BrainImg = smooth3(Img,'box',filtsizevox);

filtsizemm = 10;
filtsizevox = filtsizemm*ones(1,3).*1./XYZScale;
idx = mod(filtsizevox,2)<1;
filtsizevox = floor(filtsizevox);
filtsizevox(idx) = filtsizevox(idx)+1;

se = strel('cuboid',filtsizevox);
Img = imdilate(Img,se);
Img = imerode(Img,se);

ProjImg = smooth3(Img,'box',filtsizevox);

% DSFactor = round(mean(size(Img))/100)-2; DSFactor(DSFactor<1) = 1; DSFactor(DSFactor>3) = 3; %downsample factor (larger number for larger image)
DSFactor = 3; %going to fix this value at 3 since the mm/vox resolution is also fixed at 0.4

[BrainImgX,BrainImgY,BrainImgZ,BrainImg] = reducevolume(BrainImg,DSFactor); %downsample based on image size
BrainSurf = isosurface(BrainImgX,BrainImgY,BrainImgZ,BrainImg,0.3); %increase this value to erode image (default 0.3)
BrainSurf.vertices(:,1:2) = BrainSurf.vertices(:,[2,1]); % x and y are swapped in spm world
BrainSurf.vertices = BrainSurf.vertices*ImgInfo.mat(1:3,1:3)'+repmat(ImgInfo.mat(1:3,4)',size(BrainSurf.vertices,1),1); % rotate/scale/translate to world space

[ProjImgX,ProjImgY,ProjImgZ,ProjImg] = reducevolume(ProjImg,DSFactor); %downsample based on image size
ProjSurf = isosurface(ProjImgX,ProjImgY,ProjImgZ,ProjImg,0.3); %increase this value to erode image (default 0.3)
ProjSurf.vertices(:,1:2) = ProjSurf.vertices(:,[2,1]); % x and y are swapped in spm world
ProjSurf.vertices = ProjSurf.vertices*ImgInfo.mat(1:3,1:3)'+repmat(ImgInfo.mat(1:3,4)',size(ProjSurf.vertices,1),1); % rotate/scale/translate to world space

ProjSurf = LeG_removeSubSurf(ProjSurf); %removes any smaller enclosed surfaces within the main projection surface (i.e. ventricle remnants, etc.)
