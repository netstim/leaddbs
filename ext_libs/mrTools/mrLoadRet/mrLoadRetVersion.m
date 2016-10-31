function [ver, expectedMatlabVersion,expectedToolboxNames] = mrLoadRetVersion
%     $Id$	

ver = 4.7;

% Change this after testing Matlab upgrades
expectedMatlabVersion = [7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 7.10 7.11 7.12 7.13 7.14 8.0 8.1 8.2];

% expected toolbox
if verLessThan('matlab','8.5')
  expectedToolboxNames = {'Statistics Toolbox','Image Processing Toolbox','Optimization Toolbox'};
else
  expectedToolboxNames = {'Statistics and Machine Learning Toolbox','Image Processing Toolbox','Optimization Toolbox'};
end


