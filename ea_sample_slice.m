function [slice,boundbox]=ea_sample_slice(vol,tracor,wsize,coords,el)
interpfactor=2;

allc=[];
for side=1:length(coords)
    allc=[allc;coords{side}];
end
coords=allc;

switch tracor
    case 'tra'
        boundbox{1}=coords(el,1)-wsize:1/interpfactor:coords(el,1)+wsize;
        boundbox{2}=coords(el,2)-wsize:1/interpfactor:coords(el,2)+wsize;
        boundbox{3}=repmat(coords(el,3),1,length(boundbox{1}));
        [cmesh.X,cmesh.Y]=meshgrid(boundbox{1},boundbox{2});
        cmesh.Z=repmat(boundbox{3}(1),numel(cmesh.X),1);
        ima=spm_sample_vol(vol,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),3);
        slice=reshape(ima,length(boundbox{1}),length(boundbox{1}));
        slice=flipud(slice);
    case 'cor'
        boundbox{1}=coords(el,1)-wsize:1/interpfactor:coords(el,1)+wsize;
        boundbox{2}=repmat(coords(el,2),1,length(boundbox{1}));
        boundbox{3}=coords(el,3)-wsize:1/interpfactor:coords(el,3)+wsize;
        [cmesh.X,cmesh.Z]=meshgrid(boundbox{1},boundbox{3});
        cmesh.Y=repmat(boundbox{2}(1),numel(cmesh.X),1);
        ima=spm_sample_vol(vol,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),3);
        
        slice=reshape(ima,length(boundbox{1}),length(boundbox{1}));
        slice=flipud(slice);  
    case 'sag'
        boundbox{2}=coords(el,2)-wsize:1/interpfactor:coords(el,2)+wsize;
        boundbox{3}=coords(el,3)-wsize:1/interpfactor:coords(el,3)+wsize;
                boundbox{1}=repmat(coords(el,1),1,length(boundbox{2}));

        [cmesh.Y,cmesh.Z]=meshgrid(boundbox{2},boundbox{3});
        cmesh.X=repmat(boundbox{1}(1),numel(cmesh.Y),1);
        ima=spm_sample_vol(vol,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),3);
        
        slice=reshape(ima,length(boundbox{1}),length(boundbox{1}));
        slice=flipud(slice); 
end