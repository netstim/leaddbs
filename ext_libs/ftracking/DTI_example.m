% Example of the DTI.m Diffusion Tension Imaging (DTI) function.

% Clean everything
clear all; close all; clc

% Make a struct to store all DTI data
DTIdata=struct();

%  The test data used is from the opensource QT Fiber-Tracking, NLM insight registration & Segmentation Toolkit)

% Magnetic Gradients of data volumes 
H=[0 0 0;1 0 1;-1 0 1;0 1 1;0 1 -1; 1 1 0;-1 1 0];

%  Read the MRI (DTI) voxeldata volumes
for i=1:7,
    info = gipl_read_header(['B' num2str(i-1) '-fil.gipl']);
    DTIdata(i).VoxelData = single(gipl_read_volume(info)); 
    keyboard
    DTIdata(i).Gradient = H(i,:);
    DTIdata(i).Bvalue=1000;
end

% Constants DTI
parametersDTI=[];
parametersDTI.BackgroundTreshold=150;
parametersDTI.WhiteMatterExtractionThreshold=0.10;
parametersDTI.textdisplay=true;

% Perform DTI calculation
[ADC,FA,VectorF,DifT]=DTI(DTIdata,parametersDTI);

% Show the DiffusionTensor
figure, 
subplot(3,3,1), imshow(squeeze(DifT(:,:,round(end/2),1)),[min(DifT(:)) max(DifT(:))]); title('Dxx');
subplot(3,3,2), imshow(squeeze(DifT(:,:,round(end/2),2)),[min(DifT(:)) max(DifT(:))]); title('Dxy');
subplot(3,3,3), imshow(squeeze(DifT(:,:,round(end/2),3)),[min(DifT(:)) max(DifT(:))]); title('Dxz');
subplot(3,3,5), imshow(squeeze(DifT(:,:,round(end/2),4)),[min(DifT(:)) max(DifT(:))]); title('Dyy');
subplot(3,3,6), imshow(squeeze(DifT(:,:,round(end/2),5)),[min(DifT(:)) max(DifT(:))]); title('Dyz');
subplot(3,3,9), imshow(squeeze(DifT(:,:,round(end/2),6)),[min(DifT(:)) max(DifT(:))]); title('Dzz');

% Show the Fractional Anistropy, overlayed with the anistropy vector field.
figure, 
imshow(imresize(FA(:,:,round(end/2)),4),[]); hold on;
VectorPlotZ=squeeze(VectorF(:,:,round(end/2),1:2));
[VectorPlotX,VectorPlotY]=meshgrid(1:size(VectorPlotZ,1),1:size(VectorPlotZ,2));
quiver(VectorPlotX*4,VectorPlotY*4,VectorPlotZ(:,:,2),VectorPlotZ(:,:,1));
title('Fractional Anistropy, and Vector Field');

% Save the resulting data for the FT_test.m script.
save('FT_data','FA','VectorF');

