function [res, errStr, oArg]= dtdstruct_export(varargin)
%function [res, errStr, oArg]= dtdstruct_export(dtdStruct, comStr, [, op1[, op2[... [, opN]]]])
%
%   command:    {'enSight' | '' | 
%
%   'enSight': destDirStr, nameStr, [transMy, [dtdModStrLst]]
%           transf_My: 4x4 Matrix; M*vox -> world system
%
% Bjoern W. Kreher
% 08/05
%
% UNIX


res= []; errStr= ''; oArg= [];
maxArg= 10;

if (nargin >= 1) && dtdstruct_istype(varargin{1})
    dtdStruct= varargin{1};
else
    errStr= strcat(mfilename, ' (error): First argument have to be of the type dtdStruct');
    return
end

if (nargin >= 2) && ischar(varargin{2})
    comStr= varargin{2};
else
    errStr= strcat(mfilename, ' (error): Second argument have to be of the type string');
    return
end

argCell= cell(1, maxArg);
for i= 1:(nargin - 2)
    argCell{i}= varargin{2 + i};
end


%%%%%%%%%%%%%%%%%%%
if strcmp(comStr, 'enSight')
    if isempty(argCell{1})
        dirStr= uigetdir(pwd);
    else
        dirStr= argCell{1};
    end
    if ~exist(dirStr, 'dir')
        errStr= strcat(mfilename, ' (error): Path in dirStr does not exist');
        return    
    end

    if isempty(argCell{2}) || ischar(argCell{2})
        nameStr= argCell{2}; 
    else
        errStr= strcat(mfilename, ' (error): nameStr should be empty or a string');
        return
    end
    if isequal(size(argCell{3}), [4 4]) || isempty(argCell{3})
        transMy= argCell{3};
    else
        errStr= strcat(mfilename, ' (error): transMy should be empty or a 4x4 matrix');
        return
    end
    nameStrLst= argCell{4};
    
    [res, errStr]= local_convert2EnSight(dtdStruct, dirStr, nameStr, transMy, nameStrLst);
    
else
    errStr= strcat(mfilename, ' (error): Command ''', comStr, ''' is not supported yet');
    return    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [res, errStr]= local_convert2EnSight(dtdStruct, dirStr, nameStr, transMy, nameLst)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_convert2EnSight(dtdStruct, dirStr, identStr, transMy, nameLst)

res= {}; errStr= '';

if isempty(identStr)
    identStr= 'dtdStruct';
end

if isempty(nameLst);
    [magLst, errStr, vecLst]= dtdstruct_query(dtdStruct, 'getComStr');
    nameLst= cell(length(magLst) + length(vecLst), 1);
    nameLst(1:length(magLst))= magLst(:);
    nameLst(length(magLst)+ 1:end)= vecLst(:);
end
nameMagDataLst= {}; nameVecDataLst= {}; 

% create geo file
geo_fName= strcat(identStr, '.geo');
msgStr= strcat(mfilename, '::local_convert2EnSight (msg): exporting geometic infomation to EnSight');
disp(msgStr);
mrStruct= dtdstruct_query(dtdStruct, strcat('get', nameLst{1}));
[tmp, errStr]= mrstruct_ensight(mrStruct, 'geoFile', fullfile(dirStr, geo_fName), 1, [], 'create');
if isempty(tmp)
    return
end

% create data files
dataNo= length(nameLst);
for i= 1:dataNo
    msgStr= strcat(mfilename, '::local_convert2EnSight (msg): exporting ''', nameLst{i}, ''' to EnSight');
    disp(msgStr);
    mrStruct= dtdstruct_query(dtdStruct, strcat('get', nameLst{i}));
    if mrstruct_query(mrStruct, 'dimensions') == 3  %Magnituden map
        nameMagDataLst{end + 1, 1}= sprintf('%s_%s', identStr, nameLst{i});
        nameStr= sprintf('%s_%03d.mag', nameMagDataLst{end}, 0);
        commentStr= sprintf('Magnitude of data %s', nameLst{i});
        [tmp, errStr]= mrstruct_ensight(local_speed_work_around(mrStruct), 'dataVolume', ...
            fullfile(dirStr, nameStr), 1, commentStr);        
    elseif mrstruct_query(mrStruct, 'dimensions') == 4  % Vector field
        nameVecDataLst{end + 1, 1}= sprintf('%s_%s', identStr, nameLst{i});
        nameStr= sprintf('%s_%03d.vec', nameVecDataLst{end}, 0);    
        commentStr= sprintf('Vectorfield of data %s', nameLst{i});
        [tmp, errStr]= mrstruct_ensight(local_speed_work_around(mrStruct), 'dataVector', ...
            fullfile(dirStr, nameStr), 1, commentStr, [], 'L2');
    else
        errStr= strcat(mfilename, '::local_convert2EnSight (error): entry ''', nameLst{i}, ''' is not supported yet');
        return            
    end        
    
    if isempty(tmp)
        return
    end
end

% create case file
msgStr= strcat(mfilename, '::local_convert2EnSight (msg): Creating case file');
disp(msgStr);

casePathName  = fullfile(dirStr, strcat(identStr, '.case'));

%         
% open text file an write data
fidCase  = fopen(casePathName, 'wt');
fprintf(fidCase,'FORMAT\n');
fprintf(fidCase,'type:	 ensight gold\n');
fprintf(fidCase,'GEOMETRY\n');
fprintf(fidCase, 'model:\t %s\n', geo_fName);
fprintf(fidCase,'VARIABLE\n');
% scalar
for i= 1:length(nameMagDataLst)
    fprintf(fidCase, 'scalar per node:\t %s\t %s_***.mag\n', nameMagDataLst{i}, nameMagDataLst{i});
end
% vector
for i= 1:length(nameVecDataLst)
    fprintf(fidCase, 'vector per node:\t %s\t %s_***.vec\n', nameVecDataLst{i}, nameVecDataLst{i});
end


fprintf(fidCase,'TIME\n');
fprintf(fidCase,'time set:\t 1\n');
fprintf(fidCase,'number of steps:\t 1\n');
fprintf(fidCase,'filename start number:\t 0\n');
fprintf(fidCase,'filename increment:\t 1\n');
fprintf(fidCase,'time values:\n');
fprintf(fidCase,'0.0\n');
fclose(fidCase); % close file
% end generate case file


function [mrStrc, errStr]= local_speed_work_around(mrStrc)
errStr= '';

dimNo= mrstruct_query(mrStrc, 'dimensions');
sizeAy= mrstruct_query(mrStrc, 'sizeAy') ;
sizeAy= [sizeAy ones(1, 10)];
if dimNo == 3

elseif (dimNo == 4) && (sizeAy(4) == 3)
    sigAy= sign(mrStrc.dataAy(:, :, :, 3));
    mrStrc.dataAy(:, :, :, 1)= sigAy.*mrStrc.dataAy(:, :, :, 1);
    mrStrc.dataAy(:, :, :, 2)= sigAy.*mrStrc.dataAy(:, :, :, 2);
    mrStrc.dataAy(:, :, :, 3)= sigAy.*mrStrc.dataAy(:, :, :, 3);
else
    errStr= 'unrecognized simension';
end
    