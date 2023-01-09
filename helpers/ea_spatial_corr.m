function [rho,pval] = ea_spatial_corr(ROI1,ROI2,method)
% This code calculates spatial correlation between two binary / non-binary
% ROIS
ROI_1_tit=char(ROI1); % Path of the 1st image (First Image is Template) 
ROI_1=ea_load_nii(ROI_1_tit);

ROI_2_tit=char(ROI2); % Path of the 2nd image (First Image is Template) 
ROI_2=ea_load_nii(ROI_2_tit); 

%% Reslicing 2nd Image to Template Space (First Image)
ea_reslice_nii(ROI_2_tit,ROI_2_tit,ROI_1.voxsize);
ROI_2=ea_load_nii(ROI_2_tit);
%% Checking If NaN Exists (Transforming 0's in Binary Volumes to NaNs)
if isempty(find(isnan(ROI_1.img))) % If NaN 
    mini=min(ROI_1.img(:));
    if mini==0
        ROI_1.img(find(ROI_1.img==0))=NaN;
    end 
end 

if isempty(find(isnan(ROI_2.img))) % If NaN 
    mini=min(ROI_2.img(:));
    if mini==0
        ROI_2.img(find(ROI_2.img==0))=NaN;
    end 
end 
%% Transforming from Voxel to World Space 
nii_2=ROI_1;
mat_2=nii_2.mat;

hand=cell(2,1);
hand{1,1}=ROI_1_tit;
hand{2,1}=ROI_2_tit;

pcavox_cell=cell(2,1);
for kk=1:2
nii=ea_load_nii(hand{kk,1}); 
max_val=max(nii.img(:));
min_val=min(nii.img(:));
[X,Y,Z]=ind2sub(size(nii.img),find(nii.img(:)>min_val));
XYZvox=[X,Y,Z,ones(size(X,1),1)];
XYZmm=nii.mat*XYZvox';
%% Transforming Back from World Space to Voxel Space 
pcavox = zeros(nii_2.dim);
transf=inv(mat_2);
for i = 1:size(XYZmm,2)
    coors=transf*XYZmm(:,i);
    coors=floor(coors(1:3));
    pcavox(coors(1),coors(2),coors(3))=1;  
end 
pcavox_cell{kk,1}=pcavox;
end 
%% Calculating Spatial Correlation
vec1=pcavox_cell{1,1};
vec2=pcavox_cell{2,1};
[rho,pval]=corr(vec1(:),vec2(:),'Type',method,'Rows','complete');
disp("R-Value is "+string(rho)+", and p-value is "+string(pval)+". Method used is "+string(method))

end 