function [rho,pval] = ea_spatial_corr(ROI1_fname,ROI2_fname,method)
% This code calculates spatial correlation between two binary / non-binary
% ROIS
ROI_1_tit=char(ROI1_fname); % Path of the 1st image (First Image is Template) 
ROI_1=ea_load_nii(ROI_1_tit);
%r_ROI_1_tit=char("re_"+ROI_1_tit);

ROI_2_tit=char(ROI2_fname); % Path of the 2nd image (First Image is Template) 
r_ROI_2_tit=char("re_"+ROI_2_tit);

%% Reslicing 2nd Image to Template Space (First Image)
ea_reslice_nii(ROI_2_tit,r_ROI_2_tit,ROI_1.voxsize);
ROI_2=ea_load_nii(r_ROI_2_tit);
%% Transforming from Voxel to World Space 
%max_val=max(nii.img(:));
%min_val=min(nii.img(:));

max_val=max(ROI_2.img(:));
min_val=min(ROI_2.img(:));
[X,Y,Z]=ind2sub(size(ROI_2.img),find(ROI_2.img(:)>=min_val));
%[X,Y,Z]=ind2sub(size(ROI_2.img),find(~isnan(ROI_2.img(:))));
XYZvox=[X,Y,Z,ones(size(X,1),1)];
XYZmm=ROI_2.mat*XYZvox';
%% Transforming Back from World Space to Voxel Space 
mat_2=ROI_1.mat;
pcavox = NaN(ROI_1.dim);
transf=inv(mat_2);
for i = 1:size(XYZmm,2)
    coors=transf*XYZmm(:,i);
    coors=floor(coors(1:3));
    pcavox(coors(1),coors(2),coors(3))=ROI_2.img(X(i),Y(i),Z(i));  
end 
%% Calculating Spatial Correlation
vec1=ROI_1.img(:);
vec2=pcavox(:);

idxi_rmv=find(isnan(vec2));
vec1(idxi_rmv)=[];
vec2(idxi_rmv)=[];

[rho,pval]=corr(vec1,vec2,'Type',method);
disp("R-Value is "+string(rho)+", and p-value is "+string(pval)+". Method used is "+string(method))
end 