%% readLabels - read ITK-Snap label file 
%
% Florian Bernard, Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
   
function [labelNames, rgbColors, labelId]=readLabels(labelFile)
%    if ( Const.fileOperationsDebug )
        disp(['Reading label file: ' labelFile]);
   % end
   
    if(exist(labelFile, 'file'))
        labelFID=fopen(labelFile, 'r');
        textscan(labelFID, '%s', 14, 'delimiter','\n');
        labelRaw = textscan(labelFID, '%u %d %d %d %u %u %u %q', 'delimiter','\n');
        fclose(labelFID);
    else
        disp('Label File could not be found. Using static labels instead!');
        
       
       labelFID=['    0     0    0    0        0  0  0    "Clear Label"+'...
                 '    1     5  255    0        1  1  1    "SNr+STN_L"+' ...
                 '    2     0    0  255        1  1  1    "SNr+STN_R"+' ...
                 '    3   255   21    0        1  1  1    "NR_L"+' ...
                 '    4     0  255  255        1  1  1    "NR_R"+' ...
                 '    5   255  255    0        1  1  1    "AC"+' ...
                 '    6   255    0  255        1  1  1    "SNr_L"+' ...
                 '    7   255  239  213        1  1  1    "SNr_R"+' ...
                 '    8     0    0  205        1  1  1    "STN_L"+' ...
                 '    9   205  133   63        1  1  1    "STN_R"+' ...
                 '    10   210  180  140        1  1  1    "SNr+STN_LR"' ...
                 '    11   205  173    0        1  1  1    "MRI_Electrode_Artifact_L"' ...
                 '    12     0    0  128        1  1  1    "MRI_Electrode_Artifact_R"' ...
                 '    13   255   19    0        1  1  1    "Arterial System"']; 
                 % read labels 
        labelRaw = textscan(labelFID, '%u %d %d %d %u %u %u %q', 'delimiter','+');
	end
	labelId = labelRaw{1};
    labelNames=labelRaw{end};
    rgbColors=double([labelRaw{2} labelRaw{3} labelRaw{4}]) ./ 255;
    % attention: labelNames{1} corresponds to voxels with intensity 0
end