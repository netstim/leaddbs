%% SETUP_PACER - add all PaCER subdirectories to your MATLAB path
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2017
% mail@andreashusch.de, husch.andreas@chl.lu

%% RUN this file once (F5) to automatically make sure all PaCER files (including subdirectories) are in your MATLAB path! 
pacerDir = fileparts(mfilename('fullpath')); 
addpath(genpath(pacerDir));
savepath();