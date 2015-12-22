function [res, errStr]= morph_data(mask, kernelAy, kernelPos, OPStr, noVoxInKernel, forceKPos)
%function [res, errStr]= morph_data(mask, kernelAy, kernelPos, OPStr, noVoxInKernel, forceKPos)
%
% Funktion to apply any binary morphology kernel an binary images or 3D datasets.
%
% Input: 
%       -mask: mrStruct and schould contain a 2D or 3D data. All values ~=
%              zeros are set to one in mask
%       -kernelAy: Must have the same dimension as mask. The size of
%                  dimensions a free. 
%       -kernelPos: Describes the correspondence between mask and kernelAy
%       -OPStr: {'<='; '>='; '<'; '>'; '=='; '~='} it describes the
%               condition if a voxel will be set or not
%       -noVoxelInKernel: the number on the right side of OPStr
%       -forceKPos: {'yes'; 'no'} if set 'yes', the voxel X in res can only
%                   be set if it is also set in mask
%
% Example for 2D Erosion in a 3D-dataset:
%   kernelAy= zeros(3, 3, 3);
%   kernelAy(:, 2, 2)= 1; kernelAy(2, :, 2)= 1; kernelAy(2, 2, :)= 1;
%   erg= morph_data(mask, kernelAy, [2 2 2], '==', 5, 'no');
%
% Example for 2D dilataion in a 3D-dataset:
%   erg= morph_data(mask, kernelAy, [2 2 2], '>=', 1, 'no');
%
%
% Bjoern W. Kreher
%
% Medical Physics
% Dept. of  Diagnostic Radiology 
% University Hospital Freiburg
%
%
% 8/04
%
% UNIX

% input
% %rois= open('/home/kreher/Projects/multi_paper/data/Rois_and_visits_R2.mat');
% mask= rois.ROIstruct.f_test.dataAy;
% kernelAy= zeros(3, 3, 3);
% kernelAy(:, 2, 2)= 1; kernelAy(2, :, 2)= 1;  %kernelAy(2, 2, :)= 1; 
% kernelPos= [2, 2, 2];


res= []; errStr= [];

if ~exist('forceKPos') || isempty(forceKPos)
    forceKPos= 'No';
end


% prepare data
dimNo= mrstruct_query(mask, 'dimensions');

mapSizeAy= mrstruct_query(mask, 'sizeAy');
kernelSizeAy= size(kernelAy);
koiIdx= find(kernelAy ~= 0);
koiVc= reshape_index(koiIdx, kernelSizeAy) - 1;
koiNo= size(koiVc, 1);
kPosIdx= prod(double(ones(size(koiVc, 1), 1)*kernelPos == (koiVc + 1)), 2);

% Fehler abfangen
if (dimNo > 3) || (dimNo < 2)
    errStr= 'morph_image: Only 2D or 3D datasets are suported';
    return;
elseif dimNo ~= length(kernelSizeAy)
    errStr= 'morph_image: The dimensions of mask and kernelAy have to be equal';
    return;
elseif dimNo ~= length(kernelPos)
    errStr= 'morph_image: The dimensions of mask has to be equal to kernelPos';
    return;
end

% Generate operational data
opData= zeros([mapSizeAy + kernelSizeAy - 1 koiNo]);

if dimNo == 2
    for i= 1:koiNo
        opData(koiVc(i, 1) + (1:mapSizeAy(1)), koiVc(i, 2) + (1:mapSizeAy(2)), i)= ...
            double(mask.dataAy) * kernelAy(koiIdx(i));
    end
elseif dimNo == 3
    for i= 1:koiNo
        opData(koiVc(i, 1) + (1:mapSizeAy(1)), koiVc(i, 2) + (1:mapSizeAy(2)), koiVc(i, 3) + (1:mapSizeAy(3)), i)= ...
            double(mask.dataAy) * kernelAy(koiIdx(i));
    end
    
else
    errStr= 'morph_image: Only 2D or 3D datasets are suported';
    return;
end

% unset pixels/voxels if the center (kernelPos) is unset
if strcmpi(forceKPos, 'yes')
    opDataSizeAy= size(opData);
    opData= reshape(opData, [prod(opDataSizeAy(1:(end - 1))) koiNo]);
    idx= find(opData(:, kPosIdx) == 0);
    for i= 1:koiNo
        if strcmp(OPStr, '~=') || strcmp(OPStr, '!=')
            opData(idx, i)= 0;          %% da stimmt noch was nicht!!!
        else
            opData(idx, i)= NaN;
        end
    end
    opData= reshape(opData, opDataSizeAy);
elseif ~strcmpi(forceKPos, 'no')
    errStr= 'morph_image: forceKPos has to be ''yes'' or ''no''';
    return;
end    

% boolsche operation
if strcmp(OPStr, '=<') || strcmp(OPStr, '<=')    
    tmp= double(sum(double(opData), dimNo + 1) <= noVoxInKernel);
elseif strcmp(OPStr, '=>') || strcmp(OPStr, '>=')
    tmp= double(sum(double(opData), dimNo + 1) >= noVoxInKernel);
elseif strcmp(OPStr, '>')
    tmp= double(sum(double(opData), dimNo + 1) > noVoxInKernel);
elseif strcmp(OPStr, '<')
    tmp= double(sum(double(opData), dimNo + 1) < noVoxInKernel);
elseif strcmp(OPStr, '==') || strcmp(OPStr, '=')
    tmp= double(sum(double(opData), dimNo + 1) == noVoxInKernel);
elseif strcmp(OPStr, '~=') || strcmp(OPStr, '!=')
    tmp= double(sum(double(opData), dimNo + 1) ~= noVoxInKernel);
else
    errStr= 'morph_image: Undefined binary operator';
    return
end


% extract result
res= mask;
if dimNo == 2
    res.dataAy= tmp(kernelPos(1) - 1 + (1:mapSizeAy(1)), kernelPos(2) - 1 + (1:mapSizeAy(2)));
elseif dimNo == 3
    res.dataAy= tmp(kernelPos(1) - 1 + (1:mapSizeAy(1)), kernelPos(2) - 1 + (1:mapSizeAy(2)), kernelPos(3) - 1 + (1:mapSizeAy(3)));
else
    errStr= 'morph_image: Only 2D or 3D datasets are suported';
    res= [];
    return;
end

