%MRSTRUCT_CHECKIN checkin data into mrStruct.
%		[mrStruct, errorFlag] = mrstruct_checkin(mrStruct,dataAy,reshapeVc)
%
%		Function assigns a data array to a MR-Structure (mrStruct.dataAy=dataAy).
%		The modified mrStruct is returned. On errror, the mrStruct will
%		be returned without modification. The error status  is returned in the
%		second return value.
%		The function checks for compatibility. The dimensionality of the data set
%		has to be the same as the dimensionality of the mrStruct data type. If 
%		the dimensionality differes +-1, the mrStruct data type will be modified.
%		Function checks if mrStruct is a valid mrStruct.
%
%		Examples:
%		[mrStruct, errorFlag] = mrstruct_checkin(mrStruct,[]);		-> errorFlag=0
%
%		Harald Fischer
%		1/00
%
%		PC
% ----------------------- development -----------------------------------
%
% Authors:
%   Michael von Mengershausen (MVM)
%
% Projects:
%   MVM: renaming of matlab_new 'isvector' --> 'mn_isvector' (20060110)
%
%   FT: added dimenionality dim10 and dim11 in the years (2009 to 2011)
%
%   FT: added dimenionality dim12 (20140120)
%

function [mrStruct, errorFlag] = mrstruct_checkin(mrStruct,dataAy,reshapeVc)



%%%%% init and error check
errorFlag = 1;
if nargin==2 || nargin==3,
	if ~mrstruct_istype(mrStruct) || ~(isnumeric(dataAy) || islogical(dataAy)),
    	warning('input argument has wrong type');
        return;
    end;
else 
	warning('incorrect number of input arguments');
end;

if nargin==3,
	if ~mn_isvector(reshapeVc) || prod(reshapeVc)~=numel(dataAy(:)),
        warning('incorrect size of reshapeVc');
        return;
    end;
else
	reshapeVc = [];
end;
%%%%% End of: init and error check



%%%%% check for valid struct and dataAy
isValid = mrstruct_istype(mrStruct);
if ~isValid,
	warning('not a valid mrStruct');
	return;
end;
if ndims(dataAy)==0,
	warning('not a valid dataAy');
	return;
end;
%%%%% End of: check for valid struct and dataAy




%%%%% check in dataAy
%%%%% Note: more than one type of dataset is possible for a given dimensionality.
%%%%% Therefore, one type for each dimensionality is choosen below.
noOfDims   = mrstruct_query(mrStruct,'dimensions');
dataAyDim  = ndims(dataAy);
if mn_isvector(dataAy),
	dataAyDim = 1;
end;


%%% convert eventually
if noOfDims~=dataAyDim && noOfDims~=dataAyDim+1,
	warning('mrStruct assumes different dimensionality of dataAy - trying to convert');
    if dataAyDim==1,
    	typeStr = 'spectrum';
      
    elseif dataAyDim==2,
    	typeStr = 'image';
      
    elseif dataAyDim==3,
        typeStr = 'volume';
            
    elseif dataAyDim==4,
        typeStr = 'series3D';
      
	elseif dataAyDim==5,
    	typeStr = 'series3DEchos';
    
	elseif dataAyDim==6,
    	typeStr = '';
        
	elseif dataAyDim==7,
    	typeStr = '';

	elseif dataAyDim==8,
    	typeStr = '';
        
	elseif dataAyDim==9,
    	typeStr = '';
        
	elseif dataAyDim==10,
    	typeStr = '';
        
	elseif dataAyDim==11,
    	typeStr = 'DiffusionEchos2D';
	elseif dataAyDim==12, % FT20140120
    	typeStr = '';
% 	elseif dataAyDim==13, % FT20140212
%     	typeStr = '';
    else
    	warning('cannot checkin dataAy - dimensionality too large');
        return;
    end;
    
	tmpMrStruct   = mrstruct_init(typeStr);
    mrStruct.dim1 = tmpMrStruct.dim1;
    mrStruct.dim2 = tmpMrStruct.dim2;
    mrStruct.dim3 = tmpMrStruct.dim3;
    mrStruct.dim4 = tmpMrStruct.dim4;
    mrStruct.dim5 = tmpMrStruct.dim5;
    mrStruct.dim6 = tmpMrStruct.dim6;
    mrStruct.dim7 = tmpMrStruct.dim7;
    mrStruct.dim8 = tmpMrStruct.dim8;
    mrStruct.dim9 = tmpMrStruct.dim9;
    mrStruct.dim10 = tmpMrStruct.dim10;
    mrStruct.dim11 = tmpMrStruct.dim11;
    mrStruct.dim12 = tmpMrStruct.dim12; % FT20140120
%     mrStruct.dim13 = tmpMrStruct.dim13; % FT20140120
end;
%%% End of: convert eventually


%%% assign data
if isempty(reshapeVc),
	mrStruct.dataAy = dataAy;
else
    mrStruct.dataAy = reshape(dataAy,reshapeVc);
end;
errorFlag = 0;
%%% End of: assign data



%%%%% End of: check in dataAy
% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics
%%%%%