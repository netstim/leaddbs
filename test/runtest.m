% Runtest.m - add all LeadDBS subdirectories to your MATLAB path
%
% .. AUTHOR:
%       - Marx Loic, June 2019

% Run this file once (F5) to automatically make sure all leadDBS files (including subdirectories) are in your MATLAB path.
leadDBSDir = fileparts(mfilename('fullpath')); 
addpath(genpath(leadDBSDir));
addpath(genpath([getenv('ARTENOLIS_SOFT_PATH)' '/tools/spm12'])) 
savepath();