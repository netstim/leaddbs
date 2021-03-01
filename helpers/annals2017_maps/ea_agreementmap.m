function [outputMap, mask] = ea_agreementmap(inputMaps, outputFileName, writeoutMask, mode)
% Calculate agreement R-map

if ~exist('mode','var')
   mode='sum';
end

% Load input map
for i=1:numel(inputMaps)
    inputMaps{i} = ea_load_nii(inputMaps{i});
end

% Find agreeing voxels
positiveMask = inputMaps{1}.img>0;
negativeMask = inputMaps{1}.img<0;
for i=2:numel(inputMaps)
    positiveMask = positiveMask .* inputMaps{i}.img>0;
    negativeMask = negativeMask .* inputMaps{i}.img<0;
end
mask.both = logical(positiveMask + negativeMask);
mask.pos = positiveMask;
mask.neg = negativeMask;

% Initialize output map
outputMap = inputMaps{1};
outputMap.img = nan(size(outputMap.img));

% Multiply agreeing voxels
outputMap.img(positiveMask) = inputMaps{1}.img(positiveMask);
outputMap.img(negativeMask) = inputMaps{1}.img(negativeMask);
for i=2:numel(inputMaps)
    switch mode
        case {'mult','multiplication','multiply','*'}
            outputMap.img(positiveMask) = outputMap.img(positiveMask) .* inputMaps{i}.img(positiveMask);
            outputMap.img(negativeMask) = outputMap.img(negativeMask) .* inputMaps{i}.img(negativeMask);
        case {'sum','+'}
            outputMap.img(positiveMask) = outputMap.img(positiveMask) + inputMaps{i}.img(positiveMask);
            outputMap.img(negativeMask) = outputMap.img(negativeMask) - inputMaps{i}.img(negativeMask);
    end
end

% Ensure that negatively agreeing voxels have negative values
outputMap.img(negativeMask) = -abs(outputMap.img(negativeMask));

% Optionaly save agreement map to NIfTI
if exist('outputFileName','var')
    outputMap.fname = outputFileName;
    ea_write_nii(outputMap);

    % Optionaly save agreement map mask to NIfTI
    if exist('writeoutMask','var') && writeoutMask
        outputMask = outputMap;

        % Save mask for both positively and negatively agreeing voxels
        outputMask.fname = strrep(outputFileName, '.nii', '_mask.nii');
        outputMask.img =  mask.both;
        ea_write_nii(outputMask);

        % Save mask for positively agreeing voxels
        outputMask.fname = strrep(outputFileName, '.nii', '_posmask.nii');
        outputMask.img =  mask.pos;
        ea_write_nii(outputMask);

        % Save mask for negatively agreeing voxels
        outputMask.fname = strrep(outputFileName, '.nii', '_negmask.nii');
        outputMask.img =  mask.neg;
        ea_write_nii(outputMask);
    end
end
