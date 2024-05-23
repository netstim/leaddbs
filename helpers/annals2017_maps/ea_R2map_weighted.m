function [R2,R2adj]=ea_R2map_weighted(varargin)
% creates a correlation nifti file given a set of images and a
% regressor ("R map" in Horn 2017 AoN)

% ea_R2map_weighted(fis,regressor,outputname,mask,sk,weights) % sk can be 'k','s','sk'
% for smoothing and normalization options. Mask only necessary if
% choosing 'k' option.

pthresh=0.05;
itercount=1000;
regressor=varargin{2};
output=varargin{3};
h=nan;


[X,n]=ea_genX(varargin{1:5}); % accumulates all images into image matrix X.

if nargin>5
    weights=varargin{6};
else
    weights=ones(size(regressor));
end

X(~isfinite(X)) = 0;
nz=sum(logical(X),2)>0.2*size(X,2); % at least 20% of images covered by values.
nnz=sum(nz);
[mdls]=cellfun(@fitlm,cellfun(@transpose,mat2cell(X(nz,:),ones(1,nnz)),'un',0),repmat({regressor},1,nnz)',...
    repmat({'Weights'},1,nnz)',...
    repmat({weights},1,nnz)','Uniformoutput',0);


% maybe replace loop if can be done more easily:
Nvox=length(mdls);
R2ord=zeros(Nvox,1);
R2adj=zeros(Nvox,1);
ea_dispercent(0,'Iterating voxels');
for voxel=1:Nvox
    R2ord(voxel)=mdls{voxel}.Rsquared.Ordinary;
    R2adj(voxel)=mdls{voxel}.Rsquared.Adjusted;
    ea_dispercent(voxel/Nvox);
end
ea_dispercent(1,'end');

vargs=varargin(1:5);
[pth,fn]=fileparts(vargs{3});

EXP=nan(size(X,1),1);
EXP(nz)=R2ord;
vargs{3}=fullfile(pth,[ea_stripext(fn), '_Rord.nii']);
n.fname=vargs{3};
ea_exportmap(n,EXP,vargs{1:5});

EXP=nan(size(X,1),1);
EXP(nz)=R2adj;
vargs{3}=fullfile(pth,[ea_stripext(fn), '_Radj.nii']);
n.fname=vargs{3};
ea_exportmap(n,EXP,vargs{1:5});
