function [res, errStr, fName]= maskstruct_read(fNameIn)
%
%function [res, errStr, fName]= maskstruct_read(fName)
%
%   Opens a maskStruct from the filesystem
%
% fName:     should contain the path and file name. if it is empty or undefined, a file selection dialog will appear 
% 
% return values:
%  res:  contains the opend maskStruct in Ver2. If an error occurs res is empty
%  errStr: if an error occured, errStr contains a string to identify the error
%  fName:  if all went fine, fName contains the path of the opend ftrStruct
%
%
% Bjoern W. Kreher
% 12/07
%
% UNIX

res= []; errStr= ''; fName= '';

% file name selection
if (nargin == 0) || isempty(fNameIn) || (exist(fNameIn, 'dir') == 7)
    if (nargin > 0) && (exist(fNameIn, 'dir') == 7)
        curDir= pwd;
        cd(fNameIn);
    end
    [fileStr,dirStr]=uigetfile('*.mat','Load a maskStruct');
    if (nargin > 0) && (exist(fNameIn, 'dir') == 7)
        cd(curDir);
    end

    if fileStr == 0
        fNameIn= [];
    else
        fNameIn= fullfile(dirStr, fileStr);
    end
end
fName= fNameIn;

%% open file
if isempty(fNameIn)
    res = [];
    errStr = sprintf('User pressed cancel. No ROI selected.');
    return
end

[pathstr,fname,ext] = fileparts(fNameIn);
if isempty(ext)
    fNameIn = [fNameIn,'.mat'];
end
if ~exist(fNameIn, 'file') == 2
    res= [];
    errStr= sprintf('%s(error): file not found', mfilename);
    return
end
res= open(fNameIn);    

%% work around 
if isstruct(res) && isfield(res, 'ROIstruct')
    res= res.ROIstruct;
end


if maskstruct_istype(res, 'Ver1')
    [res, errStr]= local_convertV1_to_V2(res);
elseif maskstruct_istype(res) ~= 1
    res= [];
    errStr= 'maskstruct_read: File was not from type maskStruct';
end

%% check if dimensions ok
if size(maskstruct_query(res,'sizeAy'),2) >= 4
    res= [];
    errStr= 'maskstruct_read: maskStruct has 4 dimensions and therefore is not a valid mask';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_convertV1_to_V2(mStructV1)
res= []; errStr= '';

nameCell= fieldnames(mStructV1);
datCell= struct2cell(mStructV1);

sizeAy= [];
% some of the masks may have vox and edges set, others not
vox= [];
edges= [];
% isequal is too accurate because of possible rounding errors
tol= 1e-4;
for i= 1:length(nameCell)
    if mrstruct_istype(datCell{i})
       sizNew= mrstruct_query(datCell{i}, 'sizeAy');
       if ~isempty(sizeAy)
           if ~isequal(sizeAy, sizNew)
               errStr= sprintf('%s::local_convertV1_to_V2(error): masks in ROI file are not consistent (different volume sizes)', mfilename);
               res= [];
               return;
           end
       else
           sizeAy= sizNew;
       end               
       edgesNew= mrstruct_query(datCell{i}, 'edges');
       if ~isempty(edges)
           if ~isempty(edgesNew) && any(abs(edges(:)-edgesNew(:)) > tol)
               errStr= sprintf('%s::local_convertV1_to_V2(error): masks in ROI file are not consistent (different orientation)', mfilename);
               res= [];
               return;
           end
       else
           edges= edgesNew;
       end               
       voxNew= mrstruct_query(datCell{i}, 'vox');
       if ~isempty(vox)
           if ~isempty(voxNew) && any(abs(vox(:)-voxNew(:)) > tol)
               errStr= sprintf('%s::local_convertV1_to_V2(error): masks in ROI file are not consistent (different voxel sizes)', mfilename);
               res= [];
               return;
           end
       else
           vox= voxNew;
       end               
    end
end
mrProp= datCell{1};
mrProp.edges= edges;
mrProp.vox= vox;

if isempty(sizeAy)
    errStr= sprintf('%s::local_convertV1_to_V2(error): ROI file contains no reference volume', mfilename);
    res= [];
    return;
end
    
res= maskstruct_init(mrProp);

for i= 1:length(nameCell)
    [res, errStr]= maskstruct_modify(res, 'createMask', nameCell{i});
    if ~isempty(errStr)
        return
    end
    if ~isempty(datCell{i})
       [res, errStr]= maskstruct_modify(res, 'setMR-Properties', datCell{i});
       [res, errStr]= maskstruct_modify(res, 'setMask', datCell{i}, nameCell{i}); 
    end
    if ~isempty(errStr)
        return
    end    
end

