function [slice,boundbox,boundboxmm,sampleheight]=ea_sample_slice(vol,tracor,wsize,voxmm,coords,el)

% function samples a slice from nifti image based on coordinates and the
% wsize parameter (will use coordinate and sample a square that is 2xwsize
% long in distances). wsize can be given as mm or voxel distance (defined
% by voxmm parameter being either 'mm' or 'vox').
% Define parameter vol as spm volume (see spm_vol), define tracor as either
% 'tra', 'cor' or 'sag' for sampling direction. define coords as a set of
% points and el defining the point that is being sampled.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

interpfactor=2;

if strcmp(voxmm,'mm')
    % calculate distance in millimeters (wsize) back to voxels:
    probe=[0,0,0,1;wsize,0,0,1]';
    probe=vol.mat\probe;
    wsize=abs(round(probe(1,1)-probe(1,2)));
    clear probe
end

if iscell(coords)
    coords = vertcat(coords{:});
end

if length(coords)==1 % scalar input, only a height is defined. convert to mm space.
    getfullframe=1;
else
    getfullframe=0;
end

if any(isnan(coords(el,:)))
    %set all to nan and return
    %this is because there was no electrode here (nan coordinate)
    [slice,boundbox,boundboxmm,sampleheight]=deal(nan);
    return
end

switch tracor
    case 'tra'
        if getfullframe
            boundbox{1}=linspace(1,vol.dim(1),500);
            boundbox{2}=linspace(2,vol.dim(2),500);
            boundbox{3}=linspace(coords,coords,500);
        else
            boundbox{1}=mean(coords(el,1))-wsize:1/interpfactor:mean(coords(el,1))+wsize;
            boundbox{2}=mean(coords(el,2))-wsize:1/interpfactor:mean(coords(el,2))+wsize;
            boundbox{3}=repmat(mean(coords(el,3)),1,length(boundbox{1}));
        end
        [cmesh.X,cmesh.Y]=meshgrid(boundbox{1},boundbox{2});
        cmesh.Z=repmat(boundbox{3}(1),numel(cmesh.X),1);
        sampleheight=[1,1,boundbox{3}(1),1]';
        sampleheight=vol.mat*sampleheight;
        sampleheight=sampleheight(3);

        ima=spm_sample_vol(vol,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),3);
        slice=reshape(ima,length(boundbox{1}),length(boundbox{1}));
        %slice=fliplr(slice);
    case 'cor'
        if getfullframe
            boundbox{1}=linspace(1,vol.dim(1),500);
            boundbox{2}=linspace(coords,coords,500);
            boundbox{3}=linspace(1,vol.dim(3),500);
        else
            boundbox{1}=mean(coords(el,1))-wsize:1/interpfactor:mean(coords(el,1))+wsize;
            boundbox{2}=repmat(mean(coords(el,2)),1,length(boundbox{1}));
            boundbox{3}=mean(coords(el,3))-wsize:1/interpfactor:mean(coords(el,3))+wsize;
        end
        [cmesh.X,cmesh.Z]=meshgrid(boundbox{1},boundbox{3});
        cmesh.Y=repmat(boundbox{2}(1),numel(cmesh.X),1);
        sampleheight=[1,boundbox{2}(1),1,1]';
        sampleheight=vol.mat*sampleheight;
        sampleheight=sampleheight(2);

        ima=spm_sample_vol(vol,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),3);
        slice=reshape(ima,length(boundbox{1}),length(boundbox{1}));
        %slice=fliplr(slice);
    case 'sag'
        if getfullframe
            boundbox{1}=linspace(coords,coords,500);
            boundbox{2}=linspace(1,vol.dim(2),500);
            boundbox{3}=linspace(1,vol.dim(3),500);
        else
            boundbox{2}=coords(el,2)-wsize:1/interpfactor:coords(el,2)+wsize; % needs to be that order.
            boundbox{1}=repmat(coords(el,1),1,length(boundbox{2}));
            boundbox{3}=coords(el,3)-wsize:1/interpfactor:coords(el,3)+wsize;
        end
        [cmesh.Y,cmesh.Z]=meshgrid(boundbox{2},boundbox{3});
        cmesh.X=repmat(boundbox{1}(1),numel(cmesh.Y),1);

        sampleheight=[boundbox{1}(1),1,1,1]';
        sampleheight=vol.mat*sampleheight;
        sampleheight=sampleheight(1);

        ima=spm_sample_vol(vol,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),3);
        slice=reshape(ima,length(boundbox{1}),length(boundbox{1}));
        %slice=fliplr(slice);
end

boundboxmm{1}=vol.mat*[boundbox{1};ones(3,length(boundbox{1}))];
boundboxmm{1}=boundboxmm{1}(1,:);
boundboxmm{2}=vol.mat*[ones(1,length(boundbox{1}));boundbox{2};ones(2,length(boundbox{2}))];
boundboxmm{2}=boundboxmm{2}(2,:);
boundboxmm{3}=vol.mat*[ones(2,length(boundbox{1}));boundbox{3};ones(1,length(boundbox{3}))];
boundboxmm{3}=boundboxmm{3}(3,:);
