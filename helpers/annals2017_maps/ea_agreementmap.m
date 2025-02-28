function [outputMap, mask] = ea_agreementmap(inputMaps, outputFileName, writeoutMask, mode, zz)
% Calculate agreement R-map
% mode can be 'mult' or 'sum', 'avg' or 'sum/avg0...x' which creates an average/sum map across pixels that agree in
% polarity for 0 to x cases. (avg0 for instance would just be the average
% map across cases without any constraints.

if ~exist('mode','var')
   mode='mult';
end

% Load input map
for i=1:numel(inputMaps)
    inputMaps{i} = ea_load_nii(inputMaps{i});
end

% Find agreeing voxels
positiveMask = inputMaps{1}.img>0;
negativeMask = inputMaps{1}.img<0;
AllX=[];
for i=1:numel(inputMaps)
    if i>1
        positiveMask = positiveMask .* inputMaps{i}.img>0;
        negativeMask = negativeMask .* inputMaps{i}.img<0;
    end
    AllX=[AllX,inputMaps{i}.img(:)];
end
AllX(isnan(AllX))=0;
mask.both = logical(positiveMask + negativeMask);
mask.pos = positiveMask;
mask.neg = negativeMask;

% Initialize output map
outputMap = inputMaps{1}; % use for space
outputMap.img = nan(size(outputMap.img)); % set all voxels to nan


switch mode(1:3)
    
    case 'avg'
        if length(mode)>3
            minagree=str2double(mode(4:end));
        else
            minagree=size(AllX,2); % all have to agree
        end

        outputMap.img(:)=mean(AllX,2);
        pos=sum(AllX>0,2);
        neg=sum(AllX<0,2);
        outputMap.img((pos<minagree).*outputMap.img(:)>0)=nan;
        outputMap.img((neg<minagree).*outputMap.img(:)<0)=nan;

    case 'sum'
        if length(mode)>3
            minagree=str2double(mode(4:end));
        else
            minagree=size(AllX,2); % all have to agree
        end

        outputMap.img(:)=sum(AllX,2);
        pos=sum(AllX>0,2);
        neg=sum(AllX<0,2);
        outputMap.img((pos<minagree).*outputMap.img(:)>0)=nan;
        outputMap.img((neg<minagree).*outputMap.img(:)<0)=nan;

    case 'mul'
        minagree=size(AllX,2); % all have to agree (always the case for multiplications

        outputMap.img(:)=prod(abs(AllX),2);
        pos=sum(AllX>0,2);
        neg=sum(AllX<0,2);
        outputMap.img(logical((pos<minagree).*(neg<minagree)))=nan;
        outputMap.img(logical(neg))=-outputMap.img(logical(neg)); % retain sign
end

% set all zeros to nan
outputMap.img(all(AllX==0,2))=nan;


if exist('zz','var')
    switch zz
        case 'z'
            outputMap.img(~(outputMap.img==0))=ea_nanzscore(outputMap.img(~(outputMap.img==0)));
        case 'k'
            outputMap.img(~(outputMap.img==0))=ea_normal(outputMap.img(~(outputMap.img==0)));
    end
end

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
