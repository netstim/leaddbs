%MRSTRUCT_QUERY query mrStruct.
%		ret = mrstruct_query(mrStruct,queryStr, arg1)
%
%		queryStr = {min,max,dimensions,sizeAy,size_x,size_y,size_t,
%					size_z,echos,spectral,dataType,1...9,vox} (required).
%		'dimensions': number of dimensions in mrStruct.
%       'sizeAy':     size of dataAy in mrStruct.
%		'size_x...spectral': returns an integer from 1...9 indicating
%		the position of the queried dimension. Therefor queryStr='size_t'
%		and ret=3 means that mrStruct.dim3='size_z'.
%		'dataType': returns the type of mrStruct (image,volume...)
%		'1'...'9': Characters indicating the position of the required dimension.
%					 Returns 'size_y' for queryStr=='2' in the case of images.
%       'cmpDimensions': arg1 have to be an other mrStruct.
%           returns -1 - structs are not comparable
%                    0 - you have to coregister 
%                    1 - they are comparable but different number of
%                       dimension (eg. volume vs series 3D)
%                    2 - they are equivalent
%
%
%		Function returns integer or string values depending on the
%		query string. Both input arguments are required.
%
%		Examples:
%		posIndex = mrstruct_query(myExperimentStruct,'size_x');
%		typeStr  = mrstruct_query(myExperimentStruct,'dataType');
%
%		Harald Fischer
%		1/00
%
%		PC
%
%       Authors:
%       - Harald Fischer: HF
%       - Michael Mengershausen: MvM
%       - Bjoern Kreher: bwk
%
%       Projects:
%       - .... (HF)
%       - REGISTER_TOOL (MvM)
%
%       Changes:
%       - 30.06.05 MvM: REGISTER_TOOL: added "edges"



function [ret, errStr]= mrstruct_query(mrStruct,queryStr,arg1)



%%%%% init and error check
ret = [];   errStr= '';
if (nargin==2) || (nargin==3),
   if ~mrstruct_istype(mrStruct), %~isstruct(mrStruct) |  ~ischar(queryStr)
      warning('input argument have wrong type');
      return;
   end;
else
   warning('insufficient input arguments');
   return;
end;
%%%%% End of: init and error check



%%%%% query
if strcmp(queryStr,'min'),
   if ~isempty(mrStruct.dataAy),
      ret = min(mrStruct.dataAy(:));
   end;
   
   
elseif strcmp(queryStr,'max'),
   if ~isempty(mrStruct.dataAy),
      ret = max(mrStruct.dataAy(:));
   end;
   
   
elseif strcmp(queryStr,'dimensions'),
   if     strcmp(mrStruct.dim1,'unused'), ret=0;
   elseif strcmp(mrStruct.dim2,'unused'), ret=1;
   elseif strcmp(mrStruct.dim3,'unused'), ret=2;
   elseif strcmp(mrStruct.dim4,'unused'), ret=3;
   elseif strcmp(mrStruct.dim5,'unused'), ret=4;
   elseif strcmp(mrStruct.dim6,'unused'), ret=5;
   elseif strcmp(mrStruct.dim7,'unused'), ret=6;
   elseif strcmp(mrStruct.dim8,'unused'), ret=7;
   elseif strcmp(mrStruct.dim9,'unused'), ret=8;
   elseif strcmp(mrStruct.dim10,'unused'), ret=9;
   elseif strcmp(mrStruct.dim11,'unused'), ret=10;
   elseif strcmp(mrStruct.dim12,'unused'), ret=11;
   end;
  
elseif strcmp(queryStr,'sizeAy'),
    if ~isempty(mrStruct.dataAy),
       ret = size(mrStruct.dataAy);
    end;
  
elseif strcmp(queryStr,'sizeReal'),
    if ~isempty(mrStruct.dataAy),
       dimNo= mrstruct_query(mrStruct, 'dimensions');
       ret= ones(1, dimNo);
       sizeAy= size(mrStruct.dataAy);
       ret(1:length(sizeAy))= sizeAy;
    end;
elseif strcmp(queryStr,'size_x') || strcmp(queryStr,'size_y') ||...
        strcmp(queryStr,'size_z') || strcmp(queryStr,'size_t') ||...
        strcmp(queryStr,'echos')  || strcmp(queryStr,'spectral'),
   ret = local_find_dimension(mrStruct,queryStr);
   
   
%%% return size of each dimension in dataAy
elseif strcmp(queryStr,'dataSize'),
   y = 0;
   x = 0;
   z = 0;
   e = 0;
   t = 0;
   s = 0;
   d = 0;
   f = 0;
   if ~isempty(mrStruct.dataAy),
      dimIndex = local_find_dimension(mrStruct,'size_x');
      if dimIndex>0,
         x = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'size_y');
      if dimIndex>0,
         y = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'size_z');
      if dimIndex>0,
         z = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'echos');
      if dimIndex>0,
         e = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'size_t');
      if dimIndex>0,
         t = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'spectral');
      if dimIndex>0,
         s = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'diffusion');
      if dimIndex>0,
         d = size(mrStruct.dataAy,dimIndex);
      end;
      dimIndex = local_find_dimension(mrStruct,'ice_ind4');
      if dimIndex>0,
         f = size(mrStruct.dataAy,dimIndex);
      end;
   end;
   ret = [y,x,z,e,t,s,d];
%%% End of: return size of each dimension in dataAy


elseif strcmp(queryStr,'1') || strcmp(queryStr,'2') || strcmp(queryStr,'3') || ...
      strcmp(queryStr,'4') || strcmp(queryStr,'5')  || strcmp(queryStr,'6'),
   tmpStr = sprintf('mrStruct.dim%s',queryStr);
   ret = eval(tmpStr);
   %ret = eval(['mrStruct.dim'queryStr]); % hafi: retired

   
elseif strcmp(queryStr,'dataType'), % is slow
   dataTypeStr = 'image';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'imageEchos';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'volume';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'volumeEchos';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'series2D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'series2DEchos';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'series3D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'series3DEchos';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'spectrum';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'spectrum1D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'spectrum2D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'spectrum3D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'DiffusionEchos2D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'DiffusionEchos3D';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   dataTypeStr = 'Blade';
   if mrstruct_istype(mrStruct,dataTypeStr), ret = dataTypeStr; end;
   
% fill out remaining cases here 										%%%%%%%%%%
elseif strcmp(queryStr,'vox'),
   if ~isempty(mrStruct.vox),
       ret=mrStruct.vox;
   end;

%%%%%  MvM REGISTER_TOOL END %%%%%

elseif strcmp(queryStr,'edges'),
   if ~isempty(mrStruct.edges),
       ret=mrStruct.edges;
   end;

%%%%%  MvM REGISTER_TOOL END %%%%%


elseif strcmp(queryStr, 'cmpDimensions'),        %bwk
    if exist('arg1', 'var') && mrstruct_istype(arg1),
        mrStruct2= arg1;
    else
        errStr= sprintf('%s(cmpDimensions): invalid data', mfilename);
        return
    end;
    [ret, errStr]= local_cmpDimensions(mrStruct, mrStruct2);    

else
   warning('query string not recognized');
end;
%%%%% End of: query





%%%%%%%%%%%%%%%%%%%%% local functions %%%%%%%%%%%%%%%%%%%%%%%%

function [res, errStr]= local_cmpDimensions(mrS1, mrS2)     %bwk
res= 0; errStr= '';
epsi= 1e-4;


if ~isequal(size(mrS1.edges), [4 4]) || ~isequal(size(mrS2.edges), [4 4]),
    errStr= sprintf('%s(local_cmpDimensions): the one mrStruct does not containe a valid edge entry', mfilename);
    edgeCmp = -1
else
    diff= max(abs(mrS1.edges(:) - mrS2.edges(:)));
    if diff > epsi,
        errStr= sprintf('%s(local_cmpDimensions): the two structs have different orientations', mfilename);
        edgeCmp= 0;
    else
        edgeCmp= 1;        
    end;
end;


sizeAy_1= ones(1, 6); tmp= size(mrS1.dataAy);
sizeAy_1(1:length(tmp))= tmp;
sizeAy_2= ones(1, 6); tmp= size(mrS2.dataAy);
sizeAy_2(1:length(tmp))= tmp;

if isequal(sizeAy_1, sizeAy_2),
    sizeCmp= 1;
else
    if (mrstruct_istype(mrS1, 'image') || mrstruct_istype(mrS1, 'series2D') || ...
            mrstruct_istype(mrS1, 'imageEchos') || mrstruct_istype(mrS1, 'series2DEchos')) && ...
            isequal(sizeAy_1(1:2), sizeAy_2(1:2)),
        sizeCmp= 2;
    elseif isequal(sizeAy_1(1:3), sizeAy_2(1:3)),
        sizeCmp= 2;
    else
        sizeCmp= 0;
        errStr= sprintf('%s(local_cmpDimensions): the two structs have different size', mfilename);
    end;
end;
        
if (edgeCmp == -1)&&(sizeCmp == 0),
    res= -1;
elseif (edgeCmp == 0) || ((edgeCmp == 1)&&(sizeCmp == 0)),
    res= 0;
elseif sizeCmp == 2,
    res= 1;
else
    res= 2;
end;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dimIndex = local_find_dimension(mrStruct,queryStr)

if     strcmp(mrStruct.dim1,queryStr), dimIndex = 1;
elseif strcmp(mrStruct.dim2,queryStr), dimIndex = 2;
elseif strcmp(mrStruct.dim3,queryStr), dimIndex = 3;
elseif strcmp(mrStruct.dim4,queryStr), dimIndex = 4;
elseif strcmp(mrStruct.dim5,queryStr), dimIndex = 5;
elseif strcmp(mrStruct.dim6,queryStr), dimIndex = 6;
elseif strcmp(mrStruct.dim7,queryStr), dimIndex = 7;
elseif strcmp(mrStruct.dim8,queryStr), dimIndex = 8;
elseif strcmp(mrStruct.dim9,queryStr), dimIndex = 9;
elseif strcmp(mrStruct.dim10,queryStr), dimIndex = 10;
elseif strcmp(mrStruct.dim11,queryStr), dimIndex = 11;
else   dimIndex = 0; % not found
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics
