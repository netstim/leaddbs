% Example of the FT.m Fiber Tracking function.

% Clean everything
clear all; close all; clc

% Load the DTI fractional anistropy (FA) and Fiber VectorField
% First you have to run the DTI_test.m script !
load('FT_data','FA','VectorF');

% Read the Roi, through which all fibers must go (corpus callosum)
info = gipl_read_header('corpus_callosum.gipl');
Roi = gipl_read_volume(info)>0; 

% Fiber Tracking Constants
parametersFT=[];
parametersFT.FiberLengthMax=600;
parametersFT.FiberLengthMin=6;
parametersFT.DeviationAngleMax=1;
parametersFT.Step=0.4;
parametersFT.FiberTrackingComputationTreshold=0.125;
parametersFT.Sampling=2;
parametersFT.textdisplay=true;

% Perform fiber tracking
fibers=FT(FA,VectorF,Roi,parametersFT);

% Show FA
showcs3(FA*2.5), set(gcf, 'Renderer','OpenGL'); hold on;

% Plot all the fibers
for i=1:length(fibers),
    fiber=fibers{i};
    h=plot3t(fiber(1:2:end,1),fiber(1:2:end,2),fiber(1:2:end,3),0.2,'r');
    set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
end
view(3);
camlight;
material shiny