function [res, errStr, oArg]= mrstruct_camino(varargin)
%function [res, errStr, oArg]= mrstruct_camino(mrStruct, fNameOut, extStr, xAy, yAy, zAy)
%
%
% Bjoern W. Kreher
% 04/08
%
% UNIX


res= []; errStr= ''; oArg= [];

maxDim= 5;

if nargin < 2
    errStr= sprintf('%s(error): not enough arguments', mfilename);
    return;
end

mrStruct= varargin{1}; fNameOut= varargin{2}; extStr= varargin{3};


boxCell= cell(maxDim, 1);
for i= 1:(nargin - 3)
    boxCell{i}= varargin{i + 3};
end
    

if ~mrstruct_istype(mrStruct)
    errStr= sprintf('%s(error): first argument is not type of mrStruct', mfilename);
    return;
end


%% select fileformat
if extStr(1) == 'B'
    formatStr= 'ieee-be'; 
elseif extStr(1) == 'L'
    formatStr= 'ieee-le';
else
    errStr= sprintf('%s(error): fileformat ''%s'' is not supported yet (a precend L or B is necessary)', mfilename, extStr);
    return;
end
    
if strcmp(extStr(2:end), 'char')
    percStr= 'uint8';
elseif strcmp(extStr(2:end), 'byte')
    percStr= 'int8';
elseif strcmp(extStr(2:end), 'char')
    percStr= 'uint8';
elseif strcmp(extStr(2:end), 'short')
    percStr= 'int16';
elseif strcmp(extStr(2:end), 'int')
    percStr= 'int32';
elseif strcmp(extStr(2:end), 'long')
    percStr= 'uint64';
elseif strcmp(extStr(2:end), 'float')
    percStr= 'float32';
elseif strcmp(extStr(2:end), 'double')
    percStr= 'float64';
else
    errStr= sprintf('%s(error): fileformat ''%s'' is not supported yet', mfilename, extStr);
    return;
end

%% create index for access to volume
sizeAy= mrstruct_query(mrStruct, 'sizeAy');
dimNo= length(sizeAy);

sizeAy= ones(1, maxDim);
sizeAy(1:dimNo)= mrstruct_query(mrStruct, 'sizeAy');

newSizeAy= zeros(1, maxDim);
for i= 1:maxDim
   if isempty(boxCell{i})
       boxCell{i}= 1:sizeAy(i);
   end
   newSizeAy(i)= length(boxCell{i});
end

voxNo= prod(newSizeAy);


%% prepare index
idxNew= (1:voxNo) - 1;
idxOld= ones(1, voxNo);
segNew= prod(newSizeAy);
segOld= prod(sizeAy);

for i= maxDim:-1:1
    segNew= segNew/newSizeAy(i);
    tmpIdx= floor(idxNew/segNew);
    idxNew= mod(idxNew, segNew);
    segOld= segOld/sizeAy(i);
    idxOld= idxOld + segOld*(boxCell{i}(1 + tmpIdx) - 1);
end


%% create voxel camino file
fName= strcat(fNameOut, '.', extStr);
[fHd, errStr]= fopen(fName, 'w+', formatStr);
if fHd < 0
    return;
end

count= fwrite(fHd, mrStruct.dataAy(idxOld), percStr)

fclose(fHd);






   

