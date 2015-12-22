%MRSTRUCT_INIT create new mrStruct.
%		mrStruct = mrstruct_init(typeStr,checkinDataAy,propertiesMrStruct)
%
%		Function creates a new mrStruct of the type specified by typeStr.
%		typeStr = {image,imageEchos,volume,volumeEchos,series2D,series2DEchos,
%					spectrum,spectrum1D,spectrum2D,spectrum3D} (required).
%
%		checkinDataAy is multidimensional array and optional. If supplied,
%		checkinDataAy will be assigned to the mrStruct.dataAy.
%
%		mrStruct properties like patient name, voxel size,... are copied from
%		propertiesMrStruct if supplied (optional).
%
%		Examples:
%		studyStruct = mrstruct_init('series2D');
%		studyStruct = mrstruct_init('volume',headDataAy);
%
%		Harald Fischer
%		1/00
%
%
%       Frederik Testud (years 2009 to 2011)
%       1. Added 'diffusionEchos2D' and 'diffusionEchos3D' as mrStruct types
%       2. Added mrStruct.dim10 and mrStruct.dim11 in order to have the
%       same data size as Siemens raw data
%
%       Frederik Testud 20140120
%       1. Added 'Blade' as mrStruct types
%       2. Added mrStruct.dim12 in order to have the same data size as 
%       Siemens raw data
%
%		PC
%

function mrStruct = mrstruct_init(typeStr,checkinDataAy,propertiesMrStruct)


%%%%% init and error check
mrStruct     = [];
propMrStruct = [];
dataAy       = [];
if nargin==1 || nargin==2 || nargin==3,
    
	if ~ischar(typeStr),
    	warning('input argument is not a string - using default settings');
        typeStr = 'image';
    end;
    if ~strcmp(typeStr,'image')             &&  ~strcmp(typeStr,'imageEchos')    && ...
       ~strcmp(typeStr,'volume')            &&  ~strcmp(typeStr,'volumeEchos')   && ...
       ~strcmp(typeStr,'series2D')          &&  ~strcmp(typeStr,'series2DEchos') && ...
       ~strcmp(typeStr,'series3D')          &&  ~strcmp(typeStr,'series3DEchos') && ...
       ~strcmp(typeStr,'spectrum')          &&  ~strcmp(typeStr,'spectrum1D')    && ...
       ~strcmp(typeStr,'spectrum2D')        &&  ~strcmp(typeStr,'spectrum3D')    && ...
       ~strcmp(typeStr,'diffusionEchos2D')  &&  ~strcmp(typeStr,'diffusionEchos3D') && ...
       ~strcmp(typeStr,'Blade'), % FT20140120
        
        warning('input argument not recognized - using default settings');
        typeStr = 'image';
	end;
else
	typeStr = 'image';
end;

if nargin>=2,
	if isnumeric(checkinDataAy) || islogical(checkinDataAy),
        dataAy = checkinDataAy;
    else
    	warning('second argument is not numeric and will be ignored');
    end;
end;

if nargin==3,
	if isempty(propertiesMrStruct),
	%
    elseif mrstruct_istype(propertiesMrStruct),
        propMrStruct = propertiesMrStruct;
    else
    	warning('this argument is not a MR-Structure');
    end;
end;
%%%%% End of: init and error check



%%%%% create default struct
mrStruct.dataAy     = [];     % array: data
mrStruct.memoryType = '';     % string: 'magnitude', 'phase', 'real', 'imaginary', 'complex'

mrStruct.dim1   = 'size_x'; % string: 'size_x', 'size_y', 'size_z', 'size_t', 'spectral', 'echos', 'unused'
mrStruct.dim2   = 'size_y';
mrStruct.dim3   = 'unused';
mrStruct.dim4   = 'unused';
mrStruct.dim5   = 'unused';
mrStruct.dim6   = 'unused';
mrStruct.dim7   = 'unused';
mrStruct.dim8   = 'unused';
mrStruct.dim9   = 'unused';
mrStruct.dim10   = 'unused';
mrStruct.dim11   = 'unused';
mrStruct.dim12   = 'unused'; % FT20140120
% mrStruct.dim13   = 'unused'; % FT20140212

mrStruct.vox    = [];       % vector: [width,height,depth]
mrStruct.edges  = [];       % vector: [x1,y1,z1,x2,y2,z2,x3,y3,z3]
mrStruct.orient = '';       % string: 'axial', 'coronar', 'saggital' or '' for oblique pos.
mrStruct.method = '';       % string: name of method like 'RARE', 'MUSIC' ...
mrStruct.te     = [];       % vector:
mrStruct.tr     = [];       % real:
mrStruct.ti     = [];       % vector:
mrStruct.patient  = '';     % string: name of patient

mrStruct.user   = [];        % structure: can hold any user and method specific entries
%%%%% End of: create default struct



%%%%% set data specific properties
if strcmp(typeStr,'image')
   % default
   
elseif strcmp(typeStr,'imageEchos')
   mrStruct.dim3    = 'echos';
   
   
elseif strcmp(typeStr,'volume')
   mrStruct.dim3    = 'size_z';
   
   
elseif strcmp(typeStr,'volumeEchos')
   mrStruct.dim3    = 'size_z';
   mrStruct.dim4    = 'echos';
   
   
elseif strcmp(typeStr,'series2D')
   mrStruct.dim3    = 'size_t';
   
   
elseif strcmp(typeStr,'series2DEchos'),
   mrStruct.dim3    = 'echos';
   mrStruct.dim4    = 'size_t';
   
elseif strcmp(typeStr,'series3D'),
   mrStruct.dim3    = 'size_z';
   mrStruct.dim4    = 'size_t';
   
elseif strcmp(typeStr,'series3DEchos'),
   mrStruct.dim3    = 'size_z';
   mrStruct.dim4    = 'echos';
   mrStruct.dim5    = 'size_t';
   
elseif strcmp(typeStr,'spectrum'),
   mrStruct.dim1    = 'spectral';
   mrStruct.dim2    = 'unused';
   
elseif strcmp(typeStr,'spectrum1D'),
   mrStruct.dim2    = 'spectral';
   
elseif strcmp(typeStr,'spectrum2D'),
   mrStruct.dim3    = 'spectral';
   
elseif strcmp(typeStr,'spectrum3D'),
   mrStruct.dim3    = 'size_z';
   mrStruct.dim4    = 'spectral';
   
elseif strcmp(typeStr,'DiffusionEchoes2D'),
    mrStruct.dim11    = 'diffusion';

elseif strcmp(typeStr,'DiffusionEchoes3D'),
    mrStruct.dim3    = 'size_z';
    mrStruct.dim11    = 'diffusion';
   
elseif strcmp(typeStr,'Blade'), % FT20140120
    mrStruct.dim3    = 'size_z';
    mrStruct.dim12    = 'blades';
    
end;
%%%%% End of: set data specific properties



%%%%% checkin data if desired
if ~isempty(dataAy)
   [mrStruct, errorFlag] = mrstruct_checkin(mrStruct,dataAy);
   if errorFlag
      warning('size of data set is incompatible with type of data set - no data set checkin');
   end
end
%%%%% End of: checkin data if desired



%%%%% copy properties if desired
if ~isempty(propMrStruct)
   mrStruct.ti         = propMrStruct.ti;
   mrStruct.te         = propMrStruct.te;
   mrStruct.tr         = propMrStruct.tr;
   mrStruct.vox        = propMrStruct.vox;
   mrStruct.edges      = propMrStruct.edges;
   mrStruct.method     = propMrStruct.method;
   mrStruct.orient     = propMrStruct.orient;
   mrStruct.patient    = propMrStruct.patient;
   mrStruct.user       = propMrStruct.user;
   mrStruct.memoryType = propMrStruct.memoryType;
end;
%%%%% End of: copy properties if desired

% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics