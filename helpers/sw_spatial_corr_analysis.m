%% Sweetspot Explorer Spatial Correlation Code
% Min Jae Kim (mkim@bwh.harvard.edu)

% This code performs spatial correlations between groundtruth sweet spot
% vs. iterations of sw. generated based on different methods
ROI1='groundtruth_sweetspot.nii'; %loading groundtruth sweetspot image
sw_mapper_type="Correlations"; %type of sweetspot mapepr method used
corr_method='spearman'; % correlation method used (i,e (1) spearman or (2) pearson)
num_iterations=5; %number of sweetspot iterations generated
%% 
SW_corr_vals=[]; %1st column: r-values, 2nd column:p-values
for i=1:num_iterations
ROI2=sw_mapper_type+"_"+string(i)+".nii"; 
[rho,pval]=ea_spatial_corr_backup(ROI1,ROI2,corr_method);
SW_corr_vals(i,1)=rho;
SW_coor_vals(i,2)=pval;
end 


