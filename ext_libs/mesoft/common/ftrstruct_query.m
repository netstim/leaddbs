function [res, errStr, res2]= ftrstruct_query(varargin)
%function [res, errStr, res2]= ftrstruct_query(ftrStruct, commandStr[, op1[, op2[... [, opN]]]])
%
%   General method query a infomations and data from a ftrStruct
%
%   ftrStruct: should contain a valid ftrStruct
%   command:    {'vox' | 'getVer' | 'trackParam' | 'trackTime' | 'AlgoName' | 'hMatrix' | 'fiberNames' | 'fiberStruct' | 'roiNames' |
%                'fibersNo' | 'bBox' | 'getFiberLength' | 'getEndPoints' | 'getEndPointMap' | 'getCrossPoints' | 'getVisitMap' | 'getVisitMap_slow' |
%                'getCurve_byIndex' | 'getCurveL_byIndex' | 'getCurveL' | 'getCurveRandL' | 'getCurve' | 'getCurveRand' | 
%                'inArea' | 'endPointInArea' | 'bothEndPointInArea' | 'endPointInArea_byIndex' | 'bothEndPointInArea_byIndex'}
%   op[1..end] the parameter specified by the command string
%
% return values:
%   res: The result of the method or if an error occured res is empty
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
% Description of the different commands: 
% Information about the ftrStruct: 
%  'vox'       returns the size of a voxel
%              res: [1x4] double vector containing the voxel size of the dtdStruct [voxX voxY voxZ gapZ]
%  'getVer'    returns the vewrsion code of the ftrStruct
%              res: string containing the version code
%  'trackParam'returns the contents of the trackparameter entry
%              res: struct containing some tracking parameter
%  'trackTime' returns the contents of the trackDate entry
%              res: string containing the tracking time
%  'AlgoName'  returns the version code of the used tracking algorithm
%              res string of the version code of the tracking apram
%  'hMatrix'   returns the homogene transformation to transfor voxel inices to the world coord
%              res: [4x4] double matrix
%  'fiberNames'returns the names (identifier) of all fiber subsets (fiberStruct) in th ftrStruct
%              res: cell string containing the fiber subset names
% 
% Information about a specific fiberStruct:                
%  'fiberStruct'   returns a fiberStruct (conatining idex to the tracks) by id or name
%                  op1: name or id to identify the subset
%                  res: fiberStruct of the specified subset
%  'roiNames'  returns the ROI names, which were used for the fiber selection to create the subset
%              op1: name or id to identify the subset
%              res: cell string containing the ROI names
%  'fibersNo'  determines the number of fibertracks containing in the specified subset
%              op1: name or id to identify the subset
%              res: double(int) number fib fibers in the subset
%  'bBox'      returns the bounding box of the fiber subset
%              op1: name or id to identify the subset
%              res: [2x3] double matrix containing the min and max of the fibertrack subset
% 
% Information about the specified fiber tracks:            
%  'getFiberLength'    determines the length of all fibertracks containing in the specified subset
%              op1: name or id to identify the subset
%              res: [1xN] double array containing length of the fiber tracks
%  'getEndPoints'  determines endpoint of the fibertracks containing in the specified subset
%                  op1: name or id to identify the subset
%                  res: [2*N x 3] double vector array containing coordinates of the endpoint of the fiber tracks
%  'getEndPointMap'determines endpoint of the fibertracks containing in the specified subset and creates a mask of it
%                  op1: name or id to identify the subset
%                  op2: [1 x 3] double array containing size of the mask [sizeX sizeY sizeZ]
%                  res: [sizeX sizeY sizeZ] double array containing the mask of endpoints of the fiber tracks  
%  'getVisitMap'   determines how often a voxel was visited by the fibertracks containing in the specified subset
%                  A single fibertrack visiting a voxel twice, will increment the voxel also by two
%                  op1: name or id to identify the subset
%                  op2: [1 x 3] double array containing size of the mask [sizeX sizeY sizeZ]
%                  res: [sizeX sizeY sizeZ] double array containing the frequency a voxel was visited by the fiber tracks    
%  'getVisitMap_slow'  determines how often a voxel was visited by the fibertracks containing in the specified subset.
%                      Compared to 'getVisitMap' one fiber track can increment a voxel only on time ... and ist is slower
%                      op1: name or id to identify the subset
%                      op2: [1 x 3] double array containing size of the mask [sizeX sizeY sizeZ]
%                      res: [sizeX sizeY sizeZ] double array containing the frequency a voxel was visited by the fiber tracks 
%  'getVisitMask'  creates a mask of all voxels visited by fibertracks containing in the specified subset
%                  op1: name or id to identify the subset
%                  op2: [1 x 3] double array containing size of the mask [sizeX sizeY sizeZ]
%                  op3: threshold of minimum visits
%                  res: [sizeX sizeY sizeZ] double array containing the
%                  frequency a voxel was visited by the fiber tracks 
%  'getCurveL_byIndex' restuns the fibertracks as double array specified by an index array
%                      op1: [N x 1] double array containing index to the fibertracks
%                      op2:'n' | 'N' | []:  polygons is in respect to index coordsystem (default)
%                          'v' | 'V': polygon is scaled by the voxel size
%                          'w' | 'W': polygon is expressed in respect to the worldsystem 
%                          (only supported if hMatrix is set in ftrStruct)
%                          [4 x 4] matrix: Polygon verteces are multiplied by the matrix
%                      ret: [Ny x 3] double array. All fiber tracks are written successively separated by a NaN vector
%  'getCurve_byIndex'  restuns the fibertracks as cell array specified by an index array
%                      op1: [N x 1] double array containing index to the fibertracks
%                      op2:'n' | 'N' | []:  polygons is in respect to index coordsystem (default)
%                          'v' | 'V': polygon is scaled by the voxel size
%                          'w' | 'W': polygon is expressed in respect to the worldsystem 
%                          (only supported if hMatrix is set in ftrStruct)
%                          [4 x 4] matrix: Polygon verteces are multiplied by the matrix
%                      ret: cell array containing for each fibertrack one [Nx3] double array entry
%  'getCurveL' returns fibertracks as double array of the fiber subset
%              op1: string containing name or id of the subset (fiberStruct)
%              op2:'n' | 'N' | []:  polygons is in respect to index coordsystem (default)
%                  'v' | 'V': polygon is scaled by the voxel size
%                  'w' | 'W': polygon is expressed in respect to the worldsystem 
%                  (only supported if hMatrix is set in ftrStruct)
%                  [4 x 4] matrix: Polygon verteces are multiplied by the matrix
%              ret: [Ny x 3] double array. All fiber tracks are written successively separated by a NaN vector
%  'getCurve'  returns fibertracks as cell array of the fiber subset
%              op1: string containing name or id of the subset (fiberStruct)
%              op2:'n' | 'N' | []:  polygons is in respect to index coordsystem (default)
%                  'v' | 'V': polygon is scaled by the voxel size
%                  'w' | 'W': polygon is expressed in respect to the worldsystem 
%                  (only supported if hMatrix is set in ftrStruct)
%                  [4 x 4] matrix: Polygon verteces are multiplied by the matrix
%              ret: cell array containing for each fibertrack one [Nx3] double array entry
%  'getCurveRandL' returns randomly chosen fibertracks as double array of the fiber subset
%                  op1: string containing name or id of the subset (fiberStruct)
%                  op2: double containing max number of fiber tracks which will be rondoly chosen
%                  op3:'n' | 'N' | []:  polygons is in respect to index coordsystem (default)
%                      'v' | 'V': polygon is scaled by the voxel size
%                      'w' | 'W': polygon is expressed in respect to the worldsystem 
%                      (only supported if hMatrix is set in ftrStruct)
%                      [4 x 4] matrix: Polygon verteces are multiplied by the matrix
%                  ret: [Ny x 3] double array. All fiber tracks are written successively separated by a NaN vector   
%  'getCurveRand'  returns randomly chosen fibertracks as cell array of the fiber subset
%                  op1: string containing name or id of the subset (fiberStruct)
%                  op2: double containing max number of fiber tracks which will be rondoly chosen
%                  op3:'n' | 'N' | []:  polygons is in respect to index coordsystem (default)
%                      'v' | 'V': polygon is scaled by the voxel size
%                      'w' | 'W': polygon is expressed in respect to the worldsystem 
%                      (only supported if hMatrix is set in ftrStruct)
%                      [4 x 4] matrix: Polygon verteces are multiplied by the matrix
%                  ret: cell array containing for each fibertrack one [Nx3] double array entry
%  'inArea'    determine select or remove fibers track from a sub set and return the new fiber subset
%              op1: mrStruct containing the mask for fiber selection
%              op2: string containing name or id of the subset (fiberStruct)
%              op3: double if is set to 1, all fiber will select if they touch the mask, otherwise they will be removed
%              res: fiberStruct containing the index of the selected fibertracks
%  'endPointInArea'  
%  'bothEndPointInArea'  
%  'endPointInArea_byIndex'  
%  'bothEndPointInArea_byIndex' 
% 
%
%
% Bjoern W. Kreher
% 11/02
%
% UNIX







res= [];    res2= [];    errStr= '';

if length(varargin) < 2
    errStr= 'ftrstruct_query: There have to be at least two parameters';
    return;
end

if (~isempty(varargin)) && ftrstruct_istype(varargin{1})
    ftrStruct= varargin{1};
else
    errStr= 'ftrstruct_query: First param have to be ftrStruct';
    return;
end

if ~ischar(varargin{2})
    errStr= 'ftrstruct_query: Command have to be a string';
    return;
else
    commandStr= varargin{2};
end

for i= 1:7
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end


%%%%% begin of command switching part
if strcmp(commandStr, 'vox')
    [res, errStr]= local_getVox(ftrStruct);
    return;
elseif strcmp(commandStr, 'getVer')
    [dummy, errStr, res]= ftrstruct_istype(ftrStruct);
elseif strcmp(commandStr, 'trackParam')
    [res, errStr]= local_getTrackParam(ftrStruct);
    return;
elseif strcmp(commandStr, 'trackTime')
    [res, errStr]= local_getTrackTime(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'AlgoName')
    [res, errStr]= local_getAlgoName(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'hMatrix')
    [res, errStr]= private_getHMatrix(ftrStruct, 'W');
    return;
elseif strcmp(commandStr, 'testMRDat')
    [res, errStr]= private_testMRDat(ftrStruct, param{1});
    return;
    
elseif strcmp(commandStr, 'fiberNames')
    [res, errStr]= local_getFiberNames(ftrStruct);
    return;
elseif strcmp(commandStr, 'fiberStruct')
    [res, errStr]= local_getFiberStruct(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'roiNames')
    [res, errStr]= local_getRoiNames(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'fibersNo')
    [res, errStr]= local_getFiberNo(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'bBox')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{2});
    [res, errStr]= local_getBoundBox(ftrStruct, param{1}, trans2);
    return;
elseif strcmp(commandStr, 'getFiberLength')
    [res, errStr]= local_getFiberLength(ftrStruct, param{1}, param{2});
    return;
elseif strcmp(commandStr, 'getVertNo')
    [res, errStr]= local_getVertNo(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'getEndPoints')
    [res, errStr]= local_getEndPoints(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'getEndPointMap')
    [res, errStr]= local_getEndPointMap(ftrStruct, param{1}, param{2}, param{3});
    return;
elseif strcmp(commandStr, 'getGeneralEndPointMap')
    [res, errStr]= local_getGeneralEndPointMap(ftrStruct, param{1}, param{2}, param{3});
    return;
elseif strcmp(commandStr, 'getCrossPoints')
    [res, errStr]= local_getCrossPoints(ftrStruct, param{1});
    return;
elseif strcmp(commandStr, 'getVisitMap')
    [res, errStr]= local_getVisitMap(ftrStruct, param{1}, param{2}, param{3});
    return;
elseif strcmp(commandStr, 'getVisitMap_slow')
    [res, errStr]= local_getVisitMapSlow(ftrStruct, param{1}, param{2}, param{3});
    return;
elseif strcmp(commandStr, 'getVisitMask')
    if isempty(param{3})
        param{3} = 0;
    end        
    [res, errStr]= local_getVisitMask(ftrStruct, param{1}, param{2}, param{3},param{4});
    if isempty(errStr)
        inter= mrstruct_init('volume',res);
        inter.patient = ftrStruct.patient; %inter.user = ftrStruct.user.dtdUser;
        inter.vox = ftrStruct.vox;inter.edges = ftrStruct.hMatrix;
        res = maskstruct_init(inter);
        res = maskstruct_modify(res,'createMask','visitMask_ftr');
        res = maskstruct_modify(res,'setMask',inter.dataAy,'visitMask_ftr');
    else
        return;
    end
    return;
elseif strcmp(commandStr, 'getFractionMap')
    [res, errStr]= local_getFractionMap(ftrStruct, param{1}, param{2}, param{3});
    return
elseif strcmp(commandStr, 'getCurve_byIndex')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{2});
    [res, errStr]= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, param{1}, 'cell', trans2);
    return;
elseif strcmp(commandStr, 'getCurveL_byIndex')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{2});
    [res, errStr]= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, param{1}, 'array', trans2);
    return;
    
    
    
% elseif strcmp(commandStr, 'getDirectCurveL')
%     [trans2, errStr]= private_getHMatrix(ftrStruct, param{2});
%     [res, errStr]= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, param{1}, 'array', trans2);
%     return;
elseif strcmp(commandStr, 'getCurveL')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{2});
    [res, errStr]= local_getVertex(ftrStruct, param{1}, 'array', trans2);
    return;
elseif strcmp(commandStr, 'getCurveRandL')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{3});
    [res, errStr]= local_getVertexRand(ftrStruct, param{1}, param{2}, 'array', trans2);
    return;
elseif strcmp(commandStr, 'getCurve')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{2});
    [res, errStr]= local_getVertex(ftrStruct, param{1}, 'cell', trans2);
    return;
elseif strcmp(commandStr, 'getCurveRand')
    [trans2, errStr]= private_getHMatrix(ftrStruct, param{3});
    [res, errStr]= local_getVertexRand(ftrStruct, param{1}, param{2}, 'cell', trans2);
    return;
% elseif strcmp(commandStr, 'getVoxRand')
%     [res, errStr]= local_getVertexRand(ftrStruct, param{1}, param{2}, 'voxel', 'cell');
%     return;
elseif strcmp(commandStr, 'inArea')
    if length(param)>4,
        [res, errStr]= local_getInArea(ftrStruct, param{1}, param{2}, param{3}, param{4}, param{5});
    else
        [res, errStr]= local_getInArea(ftrStruct, param{1}, param{2}, param{3}, param{4});
    end;
    return;
elseif strcmp(commandStr, 'endPointInArea')
    [res, errStr]= local_endPointInArea(ftrStruct, param{1}, param{2}, param{3}, 'one', 'name');
elseif strcmp(commandStr, 'endPointInAreaFuzzy')
    [res, errStr]= local_endPointInArea_fuzzy(ftrStruct, param{1}, param{2}, param{3} );
elseif strcmp(commandStr, 'bothEndPointInArea')
    [res, errStr]= local_endPointInArea(ftrStruct, param{1}, param{2}, param{3}, 'both', 'name');
    
elseif strcmp(commandStr, 'endPointInArea_byIndex')
    [res, errStr]= local_endPointInArea(ftrStruct, param{1}, param{2}, param{3}, 'one', 'idx');
elseif strcmp(commandStr, 'bothEndPointInArea_byIndex')
    [res, errStr]= local_endPointInArea(ftrStruct, param{1}, param{2}, param{3}, 'both', 'idx');
    
else
    errStr= sprintf('ftrstruct_query(varagrin): command ''%s'' is not implemented', commandStr);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_testMRDat(ftr, dtd)
res= -1;  errStr= '';
epsi= 1e-4;

[dtdEdges, errStr]= dtdstruct_query(dtd, 'hMatrix');
if ~isequal(size(dtdEdges), [4 4]) || ~isequal(size(ftr.hMatrix), [4 4])
    errStr= sprintf('%s::private_testMRDat(error): the dtd/mrStruct of ftrStruct does not contain a valid edges', mfilename);
    return;
end


diff= max(abs(dtdEdges(:) - ftr.hMatrix(:)));
if diff > epsi
    errStr= sprintf('%s(local_cmpDimensions): the two structs have different orientations', mfilename);
    res= 0;
else
    res= 1;
end



%%
%
%  START:
%       [res, errStr]= local_getTrackTime(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getTrackTime(ftrStruct)

errStr= '';
res= ftrStruct.trackDate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getTrackTime
%%

%%
%
%  START:
%       [res, errStr]= local_getAlgoName(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getAlgoName(ftrStruct)

errStr= '';
res= ftrStruct.algoName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getAlgoName
%%

%%
%
%  START:
%       [res, errStr]= local_getTrackParam(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getTrackParam(ftrStruct)

errStr= '';
res= ftrStruct.TrackParam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getTrackParam
%%

%%
%
%  START:
%       [res, errStr]= local_getVox(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVox(ftrStruct)

errStr= '';
res= ftrStruct.vox;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVox
%%

%%
%
%  START:
%       [res, errStr]= local_getRoiNames(ftrStruct, nameStr;)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getRoiNames(ftrStruct, nameStr)

errStr= '';
if ~isempty(nameStr)
    [id, errStr]= private_name2id(ftrStruct, nameStr);
    res= reshape(ftrStruct.fiber{id}.roiName, [numel(ftrStruct.fiber{id}.roiName) 1]);
else
    res= {};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getRoiNames
%%

%
%
%  START:
%       [res, errStr]= local_getBoundBox(ftrStruct, nameStr;)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getBoundBox(ftrStruct, nameStr, transMy)
errStr= '';

[idxF, errStr]= private_getFiberId(ftrStruct, nameStr);
vertAy= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF, 'array', transMy);
res= [min(vertAy); max(vertAy)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getBoundBox
%

%
%
%  START:
%       [res, errStr]= local_getFiberNo(ftrStruct, nameStr;)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFiberNo(ftrStruct, nameStr)

errStr= '';
res= [];
if ~isempty(nameStr)
    [id, errStr]= private_name2id(ftrStruct, nameStr);
    if ~isempty(id)
        res= length(ftrStruct.fiber{id}.curveID);
    end
else
    res= length(ftrStruct.connectCell);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFiberNo
%

%
%
%  START:
%       [res, errStr]= local_getFiberNames(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFiberNames(ftrStruct)

errStr= '';
res= {'all'};
for i= 1:length(ftrStruct.fiber)
    res{i + 1, 1}= ftrStruct.fiber{i}.name;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFiberNames
%


%
%
%  START:
%       [res, errStr]= local_getFiberStruct(ftrStruct, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFiberStruct(ftrStruct, nameStr)

errStr= '';
if ~isempty(nameStr)
    [id, errStr]= private_name2id(ftrStruct, nameStr);
    res= ftrStruct.fiber{id};
else
    res= fiberstruct_init('fiber');
    res.name= 'all';
    res.curveId= 1:length(ftrStruct.connectCell);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFiberStruct
%

    
%
%
%  START:
%       [res, errStr]= local_getVisitMap(ftrStruct, name, sizeAy, verbose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVisitMap(ftrStruct, name, sizeAy, verbose)
res= [];    errStr= '';

[idxF, errStr]= private_getFiberId(ftrStruct, name);
if isempty(idxF)
    return;
end

res= zeros(sizeAy);
segCount= zeros(length(ftrStruct.curveSegCell), 1);

for i= 1:length(idxF)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Prepare data ... %g%%', round(1000*i/length(idxF))/10));
        drawnow
    end
    segCount(abs(ftrStruct.connectCell{idxF(i)}))= segCount(abs(ftrStruct.connectCell{idxF(i)})) + 1;
end

for i= 1:length(segCount)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Generating map ... %g%%', round(1000*i/length(segCount))/10));
        drawnow
    end
    if segCount(i) > 0
        segVc= round(0.5*(ftrStruct.curveSegCell{i}(1:(end - 1), :) + ftrStruct.curveSegCell{i}(2:end, :)));
        if not(isempty(segVc)),
            for k= 1:3
                idx= segVc(:, k) < 1;
                segVc(idx, k)= 1;
                idx= segVc(:, k) > sizeAy(k);
                segVc(idx, k)= sizeAy(k);
            end
            segIdx= reshape_index_back(segVc, sizeAy);
            res(segIdx)= res(segIdx) + segCount(i);
        end;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVisitMap
%
%

%
%
%  START:
%       [res, errStr]= local_getFractionMap(ftrStruct, name, sizeAy, verbose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFractionMap(ftrStruct, name, sizeAy, verbose)
res= [];    errStr= '';

[vertCell, errStr]= local_getVertex(ftrStruct, name, 'cell', []);

if isempty(vertCell)
    return
end

res= zeros(sizeAy);
    
% fuer alle fibers ...
for i= 1:length(vertCell)
    vAy= vertCell{i};
    idx= reshape_index_back(round(0.5*(vAy(1:(end-1), :) + vAy(2:end, :))), sizeAy);
    lenAy= sqrt(sum((vAy(1:(end-1), :) - vAy(2:end, :)).^2, 2));
    for ii= 1:length(idx)
        res(idx(ii))= res(idx(ii)) + lenAy(ii); 
    end
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Process data ... %g%%', round(1000*i/length(vertCell))/10));
        drawnow;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFractionMap
%
%

%
%
%  START:
%       [res, errStr]= local_getVisitMapSlow(ftrStruct, name, sizeAy, verbose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVisitMapSlow(ftrStruct, name, sizeAy, verbose)
res= [];    errStr= '';

[idxF, errStr]= private_getFiberId(ftrStruct, name);

if isempty(idxF)
    return
end

res= zeros(sizeAy);
segCount= zeros(length(ftrStruct.curveSegCell), 1);

%%% mit fibergewichtung oder nicht
if isempty(ftrStruct.user)
    userFields= [];
else
    userFields= fieldnames(ftrStruct.user);
end
segFrak= 0;
for i= 1:length(userFields)
    if strcmp(userFields{i}, 'segFraction')
        segFrak= 1;
    end
end
    
% fuer alle fibers ...
for i= 1:length(idxF)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Process data ... %g%%', round(1000*i/length(idxF))/10));
        drawnow
    end
    seg_Ids= abs(ftrStruct.connectCell{idxF(i)});
    tmpAy= zeros(sizeAy);
    for curSeg=seg_Ids
        segVc= round(0.5*(ftrStruct.curveSegCell{curSeg}(1:(end - 1), :) + ftrStruct.curveSegCell{curSeg}(2:end, :)));
        for k= 1:3
            idx= segVc(:, k) < 1;
            segVc(idx, k)= 1;
            idx= segVc(:, k) > sizeAy(k);
            segVc(idx, k)= sizeAy(k);
        end
        segIdx= reshape_index_back(segVc, sizeAy);
        if segFrak
            tmpAy(segIdx)= max(tmpAy(segIdx), ones(size(segIdx))*ftrStruct.user.segFraction(curSeg));
        else
            tmpAy(segIdx)= 1;
        end
    end
    res= res + tmpAy;
end

%res= 100*res/length(idxF);
%res= 100*res/length(idxF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVisitMapSlow
%
%

%%
%
%
%  START:
%       [res, errStr]= local_getVisitMask(ftrStruct, name, sizeAy, verbose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVisitMask(ftrStruct, name, sizeAy, thresh, verbose)
res= [];    errStr= '';

[idxF, errStr]= private_getFiberId(ftrStruct, name);
if isempty(idxF)
    return;
end

res= zeros(sizeAy);
segCount= zeros(length(ftrStruct.curveSegCell), 1);

for i= 1:length(idxF)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Prepare data ... %g%%', round(1000*i/length(idxF))/10));
        drawnow
    end
    segCount(abs(ftrStruct.connectCell{idxF(i)}))= segCount(abs(ftrStruct.connectCell{idxF(i)})) + 1;
end

for i= 1:length(segCount)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Generating map ... %g%%', round(1000*i/length(segCount))/10));
        drawnow
    end
    if segCount(i) > 0
        segVc= round(0.5*(ftrStruct.curveSegCell{i}(1:(end - 1), :) + ftrStruct.curveSegCell{i}(2:end, :)));
        for k= 1:3
            idx= segVc(:, k) < 1;
            segVc(idx, k)= 1;
            idx= segVc(:, k) > sizeAy(k);
            segVc(idx, k)= sizeAy(k);
        end
        segIdx= reshape_index_back(segVc, sizeAy);
        res(segIdx)= res(segIdx) + segCount(i);

    end
end
res(res <= thresh) = 0;
res(res > thresh) = 1;
res = logical(res);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVisitMask

%%

%%
%
%  START:
%   [res, errStr]= local_getVertNo(ftrStruct, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVertNo(ftrStruct, name)

res= [];    errStr= '';
[idxF, errStr]= private_getFiberId(ftrStruct, name);

if isempty(idxF)
    errStr= 'local_getVertNo: unknown fiber subset';
    return;
end

fiberCell= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF, 'cell');
res= zeros(1, length(fiberCell));

for i= 1:length(fiberCell)
    res(i)= size(fiberCell{i}, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVertNo
%
%



%
%
%  START:
%   [res, errStr]= local_getFiberLength(ftrStruct, name, verbose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFiberLength(ftrStruct, name, verbose)
res= [];    errStr= '';
[idxF, errStr]= private_getFiberId(ftrStruct, name);

if isempty(idxF)
    return;
end

if ~isempty(verbose) 
    set(verbose, 'String', sprintf('getFiberLength: prepare data of fiber ''%s''', name)); 
    drawnow
end

if length(ftrStruct.vox) == 3
    vox= reshape(ftrStruct.vox, [1 3]);
elseif length(ftrStruct.vox) == 4
    vox= [ftrStruct.vox(1) ftrStruct.vox(2) ftrStruct.vox(3) + ftrStruct.vox(4)];
else
    errStr= 'No vox size is defined in ftrStruct';
    return
end

fiberCell= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF, 'cell');
res= zeros(1, length(fiberCell));

for i= 1:length(fiberCell)
    vertVc= fiberCell{i}.*(ones(size(fiberCell{i}, 1), 1)*vox);
    distAy= sqrt(sum((vertVc(1:(end - 1), :) - vertVc(2:end, :)).^2, 2));
    res(i)= sum(distAy);
    if ~isempty(verbose) && (mod(i, 500) == 0)
        set(verbose, 'String', sprintf('getFiberLength: calculation length of fiber ''%s'' (%g%% done)', name, round(100*i/length(fiberCell))));
        drawnow
    end
end

if ~isempty(verbose) 
    set(verbose, 'String', 'getFiberLength: done'); 
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFiberLength
%
%

%
%
%  START:
%       [res, errStr]= local_getEndPointMap(ftrStruct, name)
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getGeneralEndPointMap(ftrStruct, name, sizeAy, verbose)
res= [];    errStr= '';
[idxF, errStr]= private_getFiberId(ftrStruct, name);

if isempty(idxF)
    return;
end

if ~isempty(verbose)
    set(verbose, 'String','Prepare data for EndPoint-Filter ...');
    drawnow
end

res= zeros(sizeAy);
fiberCell= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF, 'cell');
for i= 1:length(fiberCell)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Generating EndPoint-Map ... %g%%', round(1000*i/length(idxF))/10));
        drawnow
    end
    epVc= round(fiberCell{i}([1 end], :));
    res(epVc(1, 1), epVc(1, 2), epVc(1, 3))= res(epVc(1, 1), epVc(1, 2), epVc(1, 3)) + 1;
    res(epVc(2, 1), epVc(2, 2), epVc(2, 3))= res(epVc(2, 1), epVc(2, 2), epVc(2, 3)) + 1;
end




%
%
%  START:
%       [res, errStr]= local_getEndPointMap(ftrStruct, name)
% Nur f?r FACT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getEndPointMap(ftrStruct, name, sizeAy, verbose)
res= [];    errStr= '';
[idxF, errStr]= private_getFiberId(ftrStruct, name);

if isempty(idxF)
    return;
end

if ~isempty(verbose)
    set(verbose, 'String','Prepare data for EndPoint-Filter ...');
    drawnow
end

epVc= zeros(2, 3);  rVert= zeros(2, 3);
res= zeros(sizeAy);
maxMy= (ones(2, 1)*sizeAy);
for i= 1:length(idxF)
    if (mod(i, 100) == 0) && ~isempty(verbose)
        set(verbose, 'String', ...
            sprintf('Generating EndPoint-Map ... %g%%', round(1000*i/length(idxF))/10));
        drawnow
    end
    segId= ftrStruct.connectCell{idxF(i)}(1);
    if segId > 0
        firstVc= ftrStruct.curveSegCell{segId}(1, :);
        rVert(1, :)= round(0.5*(firstVc + ftrStruct.curveSegCell{segId}(2, :)));
    else
        firstVc= ftrStruct.curveSegCell{-segId}(end, :);
        rVert(1, :)= round(0.5*(firstVc + ftrStruct.curveSegCell{-segId}(end - 1, :)));
    end
    
    segId= ftrStruct.connectCell{idxF(i)}(end);
    if segId > 0
        lastVc= ftrStruct.curveSegCell{segId}(end, :);
        rVert(2, :)= round(0.5*(lastVc + ftrStruct.curveSegCell{segId}(end - 1, :)));
    else
        lastVc= ftrStruct.curveSegCell{-segId}(1, :);
        rVert(2, :)= round(0.5*(lastVc + ftrStruct.curveSegCell{-segId}(2, :)));
    end

    epVc(1, :)= round(firstVc - rVert(1, :)) + rVert(1, :);
    epVc(2, :)= round(lastVc - rVert(2, :)) + rVert(2, :);

    epVc= (epVc < 1) + (epVc >= 1).*epVc;
    epVc= (epVc >= maxMy).*maxMy + (epVc < maxMy).*epVc;
    res(epVc(1, 1), epVc(1, 2), epVc(1, 3))= res(epVc(1, 1), epVc(1, 2), epVc(1, 3)) + 1;
    res(epVc(2, 1), epVc(2, 2), epVc(2, 3))= res(epVc(2, 1), epVc(2, 2), epVc(2, 3)) + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getEndPointMap
%
%

%
%
%  START:
%       [res, errStr]= local_getEndPoints(ftrStruct, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getEndPoints(ftrStruct, name)
errStr= ''; res= [];

[idxF, errStr]= private_getFiberId(ftrStruct, name);
if isempty(idxF), return; end

res= zeros(2*length(idxF), 3);
for i= 1:length(idxF)
    segId= ftrStruct.connectCell{idxF(i)}(1);
    if segId > 0
        res(2*i - 1, :)= ftrStruct.curveSegCell{segId}(1, :);
    else
        res(2*i - 1, :)= ftrStruct.curveSegCell{-segId}(end, :);
    end       
    segId= ftrStruct.connectCell{idxF(i)}(end);
    if segId > 0
        res(2*i, :)= ftrStruct.curveSegCell{segId}(end, :);
    else
        res(2*i, :)= ftrStruct.curveSegCell{-segId}(1, :);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getEndPoints
%
%

%
%  START:
%       [res, errStr]= local_getCrossPoints(ftrStruct, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getCrossPoints(ftrStruct, name, type)
errStr= ''; res= [];
[idxF, errStr]= private_getFiberId(ftrStruct, name);
if isempty(idxF), return; end

crossC= 0;
for i= 1:length(idxF)
    crossC= crossC + length(ftrStruct.connectCell{idxF(i)}) - 1;
end
res= zeros(crossC, 3);
count= 1;
for i= 1:length(idxF)
    idxS= ftrStruct.connectCell{idxF(i)};
    if length(idxS) > 1
        segId= idxS(1);
        if segId > 0
            pervVc= ftrStruct.curveSegCell{segId}(end, :);
        else
            pervVc= ftrStruct.curveSegCell{-segId}(1, :);
        end
        
        for j= 2:length(idxS)
            segId= idxS(j);
            if segId > 0
                curVc= ftrStruct.curveSegCell{segId}(1, :);
                nextVc= ftrStruct.curveSegCell{segId}(end, :);
            else
                curVc= ftrStruct.curveSegCell{-segId}(end, :);
                nextVc= ftrStruct.curveSegCell{-segId}(1, :);
            end
            res(count, :)= round(.5*(pervVc + curVc));
            prevVc= nextVc;
            count = count + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getCrossPoints
%

%
%
%  START:
%       [res, errStr]= local_getVertex(ftrStruct, ratio)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVertex(ftrStruct, name, type, transMy)

errStr= ''; res=  [];
[idxF, errStr]= private_getFiberId(ftrStruct, name);
if isempty(idxF), return; end
res= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF, type, transMy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVertex
%

%
%
%  START:
%       [res, errStr]= local_getVertexRand(ftrStruct, ratio)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVertexRand(ftrStruct, name, cNo, type, transMy)

errStr= ''; res=  [];

[idxF, errStr]= private_getFiberId(ftrStruct, name);
if isempty(idxF), return; end

ratio= cNo/length(idxF);
rndTmp= rand(1, length(idxF)) < ratio;
res= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF(rndTmp), type, transMy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVertexRand
%




%
%
%  START:
%       [res, errStr]= local_getInArea(ftrStruct, region, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fiberStruct, errStr]= local_getInArea_old(ftrStruct, region, nameRef, not, txtHandle)
if exist('txtHandle')
    verbose= 1;
else
    verbose= 0;
end

if verbose
    set(txtHandle, 'String', 'Prepare data ...');
    drawnow;
end
[idxF, errStr]= private_getFiberId(ftrStruct, nameRef);
vertices= private_getVertex(ftrStruct.connectCell, ftrStruct.curveSegCell, idxF, 'cell');

idxF_new= [];

if verbose
    set(txtHandle, 'String', 'Starting ...');
    drawnow;
end

roiVert= reshape_index(find(region.dataAy), mrstruct_query(region, 'sizeAy'));
for i= 1:length(vertices)
    vertex= round(0.5*(vertices{i}(2:end, :) + vertices{i}(1:(end-1), :)));
    if ~isempty(roiVert) && ~isempty(vertex)
        cmp= zeros(size(vertex, 1), size(roiVert, 1));
        for l= 1:3
            [a, b]= meshgrid(roiVert(:, l), vertex(:, l));
            cmp= cmp + (a == b);
        end
        if ~any(cmp(:) == 3)
            idxF_new(end + 1)= idxF(i);
        end
    end
    if verbose && (mod(i, 100) == 0)
        if not
            set(txtHandle, 'String', sprintf('%3.1f%% curves finished and %d fibers found until now', ...
                round(1000*i/length(vertices))/10, length(idxF_new)));
        else
            set(txtHandle, 'String', sprintf('%3.1f%% curves finished and %d fibers cut out until now', ...
                round(1000*i/length(vertices))/10, length(idxF_new)));
        end
        drawnow;
    end
end

if isempty(idxF_new)
    fiberStruct= [];
    errStr= 'ftrstruct_query (local_getInArea): No fibers found';
else
    fiberStruct= ftrstruct_init('fiber');
    fiberStruct.curveID= idxF_new;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVertexRand
%


%
%
%  START:
%       [res, errStr]= local_getInArea(ftrStruct, region, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fiberStruct, errStr]= local_getInArea(ftrStruct, region, nameRef, not, txtHandle, fuzzy)
if ~exist('txtHandle') || isempty(txtHandle)
    verbose= 0;
else
    verbose= 1;
end

if verbose
    set(txtHandle, 'String', 'Prepare data ...');
    drawnow;
end
[idxF, errStr, fstr]= private_getFiberId(ftrStruct, nameRef);

if verbose
    set(txtHandle, 'String', 'Starting ...');
    drawnow;
end

if ~exist('fuzzy') || isempty(fuzzy)
    isfuzzy= 0;
else
    isfuzzy= 1;
end

segFlag= uint8(zeros(length(ftrStruct.curveSegCell), 1));
idxF_new= [];
%%% vergleiche welche Segmente in Region sind
%roiVert= reshape_index(find(region.dataAy), mrstruct_query(region, 'sizeAy'));


set(txtHandle, 'String', 'checking'); drawnow;
segIdx= cat(1,ftrStruct.connectCell{idxF(:)});
vertex =  (cat(1,ftrStruct.curveSegCell{segIdx(:)}));
sz = size(region.dataAy);
vertex(vertex(:,1)<1,1) = 1;
vertex(vertex(:,2)<1,2) = 1;
vertex(vertex(:,3)<1,3) = 1;
vertex(vertex(:,1)>sz(1),1) = sz(1);
vertex(vertex(:,2)>sz(2),2) = sz(2);
vertex(vertex(:,3)>sz(3),3) = sz(3);
lens = cellfun(@(x) size(x,1), ftrStruct.curveSegCell(segIdx(:)));
cumlens = cumsum(lens);
region.dataAy = double(region.dataAy);


thres = 0;
if ~isfuzzy,
    vertex = round(vertex);        
    idx = vertex(:,1) + sz(1)*(vertex(:,2)-1) + (sz(1)*sz(2))*(vertex(:,3)-1);
    cmpall = region.dataAy(idx);        
else       
    if fuzzy.sigma > 0,
       sz = size(region.dataAy);
       [X Y Z] = ndgrid(-ceil(sz(1)/2):floor(sz(1)/2)-1,-ceil(sz(2)/2):floor(sz(2)/2)-1,-ceil(sz(3)/2):floor(sz(3)/2)-1);
       X = fftshift(X); Y = fftshift(Y); Z = fftshift(Z);
       R2 = X.^2 + Y.^2 + Z.^2;
       isog = fftn(exp(-R2/(2*fuzzy.sigma^2)));       
       isog = isog/max(isog(:));
       region.dataAy = real(ifftn(fftn(squeeze(region.dataAy)).*isog));
    end;    
    
    cmpall = interp3(region.dataAy,vertex(:,2),vertex(:,1),vertex(:,3));                
    thres = fuzzy.thres;
end;

set(txtHandle, 'String', 'assembling result'); drawnow;
cumcmpall = cumsum(cmpall);
overlap = (cumcmpall([cumlens]) - [0 ; cumcmpall(cumlens(1:end-1))]);

if isfield(fuzzy,'relative'),
    if fuzzy.relative
        overlap = overlap(:) ./(eps+lens(:));
    end;
end;

if not == 1,
     locidx = overlap>thres;
else
     locidx = overlap<=thres;
end;
idxF_new = idxF(locidx);
overlap = overlap(locidx);


    
    
%%
% 
% cnt = 1;
% for k= 1:length(idxF)
%     segIdx= ftrStruct.connectCell{idxF(k)};
%     contains= 0;
%     for i= 1:length(segIdx)
%         if segFlag(segIdx(i)) == 2
%             contains= 1;
%         elseif  segFlag(segIdx(i)) == 0
%             
%             len =  size(ftrStruct.curveSegCell{segIdx(i)},1);
%             if len > 0,
%                 cmp = cmpall(cnt:cnt+len-1);
%                 if any(cmp(:)>0)
%                     contains= 1;
%                     segFlag(segIdx(i))= 2;
%                 else
%                     segFlag(segIdx(i))= 1;
%                 end                
%             end
%             cnt = cnt +len;
%             
%         end
%     end
%     if contains == not
%         idxF_new(end + 1)= idxF(k);
%     end
%     
%     if (mod(k, round(length(idxF)/20)) == 0)
%         verbStr= sprintf('Comparing: %3.1f%% curves finished ', round(1000*k/length(idxF))/10);
%         if verbose
%             set(txtHandle, 'String', verbStr);
%             drawnow;
%         else
%             if (mod(k, 4000) == 0)
%                 disp(verbStr);
%             end
%         end
%     end
% end



% 
% 
% for k= 1:length(idxF)
%     segIdx= ftrStruct.connectCell{idxF(k)};
%     contains= 0;
%     for i= 1:length(segIdx)
%         if segFlag(segIdx(i)) == 2
%             contains= 1;
%         elseif  segFlag(segIdx(i)) == 0
% %            vertLen= size((ftrStruct.curveSegCell{segIdx(i)}), 1) + 1;
% %            vertex= zeros(vertLen, 3);
% %            vertex(2:(vertLen-1), :)= round(0.5*(ftrStruct.curveSegCell{segIdx(i)}(2:end, :) + ftrStruct.curveSegCell{segIdx(i)}(1:(end-1), :)));
% %            vertex(1, :)= round(ftrStruct.curveSegCell{segIdx(i)}(1, :) - vertex(2, :)) + vertex(2, :);
% %            vertex(end, :)= round(ftrStruct.curveSegCell{segIdx(i)}(end, :) - vertex(end - 1, :)) + vertex(end - 1, :);
%             vertex =  round(ftrStruct.curveSegCell{segIdx(i)});
%             if ~isempty(roiVert) & ~isempty(vertex)
%                 sz = size(region.dataAy);
%                 vertex(vertex(:,1)<1,1) = 1;
%                 vertex(vertex(:,2)<1,2) = 1;
%                 vertex(vertex(:,3)<1,3) = 1;
%                 vertex(vertex(:,1)>sz(1),1) = sz(1);
%                 vertex(vertex(:,2)>sz(2),2) = sz(2);
%                 vertex(vertex(:,3)>sz(3),3) = sz(3);
%                 
%                 idx = vertex(:,1) + sz(1)*(vertex(:,2)-1) + (sz(1)*sz(2))*(vertex(:,3)-1);
%                 cmp = region.dataAy(idx);
%                 if any(cmp(:)>0)
%                     contains= 1;
%                     segFlag(segIdx(i))= 2;
%                 else
%                     segFlag(segIdx(i))= 1;
%                 end
% %                cmp= zeros(size(vertex, 1), size(roiVert, 1));
% %                 for l= 1:3
% %                     [a, b]= meshgrid(roiVert(:, l), vertex(:, l));
% %                     cmp= cmp + (a == b);
% %                 end
% %                 if any(cmp(:) == 3)
% %                     contains= 1;
% %                     segFlag(segIdx(i))= 2;
% %                 else
% %                     segFlag(segIdx(i))= 1;
% %                 end
%                 
%             end
%         end
%     end
%     if contains == not
%         idxF_new(end + 1)= idxF(k);
%     end
%     
%     if (mod(k, round(length(idxF)/20)) == 0)
%         verbStr= sprintf('Comparing: %3.1f%% curves finished ', round(1000*k/length(idxF))/10);
%         if verbose
%             set(txtHandle, 'String', verbStr);
%             drawnow;
%         else
%             if (mod(k, 4000) == 0)
%                 disp(verbStr);
%             end
%         end
%     end
% end
% 
% 







if isempty(idxF_new)
    fiberStruct= [];
    errStr= 'ftrstruct_query (local_getInArea): No fibers found';
else
    fiberStruct= ftrstruct_init('fiber');
    fiberStruct.curveID= idxF_new;
    
    if isfuzzy,
        fiberStruct.user.overlap = overlap(:);
        fiberStruct.user.thres = thres;
        if ~(isempty(fstr)) 
            if isfield(fstr.user,'overlap'),
                fiberStruct.user.overlap = [overlap(:) fstr.user.overlap(locidx,:)];
                fiberStruct.user.thres = [thres ; fstr.user.thres(:)];
            end;
        end;
    end;
    
    
    
end









%
%
%  START:
%       [res, errStr]= local_getInArea(ftrStruct, region, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [fiberStruct, errStr]= local_endPointInArea_fuzzy(ftrStruct, region, nameRef,fuzzy)

fiberStruct= [];    errStr= '';

[idxF, errStr]= private_getFiberId(ftrStruct, nameRef);
if isempty(idxF)
    return
end

if ~exist('fuzzy') || isempty(fuzzy)
    isfuzzy= 0;
else
    isfuzzy= 1;
end

segFlag= uint8(zeros(length(ftrStruct.curveSegCell), 1));
idxF_new= [];
%%% vergleiche welche Segmente in Region sind
%roiVert= reshape_index(find(region.dataAy), mrstruct_query(region, 'sizeAy'));


endPoints = cellfun(@(x) [x(1,:) ; x(end,:)], ftrStruct.curveSegCell,'uniformoutput',false);

segIdx= cat(1,ftrStruct.connectCell{idxF(:)});
vertex = (cat(1,endPoints{segIdx(:)})); %(cat(1,ftrStruct.curveSegCell{segIdx(:)}));
sz = size(region.dataAy);
vertex(vertex(:,1)<1,1) = 1;
vertex(vertex(:,2)<1,2) = 1;
vertex(vertex(:,3)<1,3) = 1;
vertex(vertex(:,1)>sz(1),1) = sz(1);
vertex(vertex(:,2)>sz(2),2) = sz(2);
vertex(vertex(:,3)>sz(3),3) = sz(3);
lens = cellfun(@(x) size(x,1), endPoints(segIdx(:)));
cumlens = cumsum(lens);
region.dataAy = double(region.dataAy);


thres = 0;
if ~isfuzzy,
    vertex = round(vertex);        
    idx = vertex(:,1) + sz(1)*(vertex(:,2)-1) + (sz(1)*sz(2))*(vertex(:,3)-1);
    cmpall = region.dataAy(idx);        
else       
    if fuzzy.sigma > 0,
       sz = size(region.dataAy);
       [X Y Z] = ndgrid(-ceil(sz(1)/2):floor(sz(1)/2)-1,-ceil(sz(2)/2):floor(sz(2)/2)-1,-ceil(sz(3)/2):floor(sz(3)/2)-1);
       X = fftshift(X); Y = fftshift(Y); Z = fftshift(Z);
       R2 = X.^2 + Y.^2 + Z.^2;
       isog = fftn(exp(-R2/(2*fuzzy.sigma^2)));       
       isog = isog/max(isog(:));
       region.dataAy = real(ifftn(fftn(squeeze(region.dataAy)).*isog));
    end;    
    
    cmpall = interp3(region.dataAy,vertex(:,2),vertex(:,1),vertex(:,3));                
    thres = fuzzy.thres;
end;

%set(txtHandle, 'String', 'assembling result'); drawnow;
cumcmpall = cumsum(cmpall);
overlap = (cumcmpall([cumlens]) - [0 ; cumcmpall(cumlens(1:end-1))]);

if isfuzzy,
    if isfield(fuzzy,'relative'),
        if fuzzy.relative
            overlap = overlap(:) ./(eps+lens(:));
        end;
    end;
end;

locidx = overlap>thres;

idxF_new = idxF(locidx);
overlap = overlap(locidx);





if isempty(idxF_new)
    fiberStruct= [];
    errStr= 'ftrstruct_query (local_getInArea): No fibers found';
else
    fiberStruct= ftrstruct_init('fiber');
    fiberStruct.curveID= idxF_new;
    
    if isfuzzy,
        fiberStruct.user.overlap = overlap(:);
        fiberStruct.user.thres = thres;
    end;
    
    
    
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVertexRand
%


%
%
%  START:
%       [res, errStr]= local_getInArea(ftrStruct, region, name)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fiberStruct, errStr]= local_endPointInArea(ftrStruct, region, nameRef, txtHandle, modeStr, adrModeStr)

fiberStruct= [];    errStr= '';

if isempty(modeStr)
    modeStr= 'one';
end

if isempty(txtHandle)
    verbose= 0;
else
    verbose= 1;
end
if strcmp(modeStr, 'both')
    bothFlag= 1;
else
    bothFlag= 0;
end
if isfield(ftrStruct.user, 'endA_Vc') && isfield(ftrStruct.user, 'endB_Vc')
    endFlag= 1;
else
    endFlag= 0;
end
% if verbose
%     set(txtHandle, 'String', 'Prepare data ...');
%     drawnow;
% end

if strcmp(adrModeStr, 'name')
    [idxF, errStr]= private_getFiberId(ftrStruct, nameRef);
else
    idxF= nameRef;
end
if isempty(idxF)
    return
end
% 
% if verbose
%     set(txtHandle, 'String', 'Starting ...');
%     drawnow;
% end

roiVert= reshape_index(find(region.dataAy), mrstruct_query(region, 'sizeAy'));

if isempty(roiVert)
    errStr= 'ftrstruct_query (local_endPointInArea): roi is empty';
    return
end

tmpVc= zeros(2, 3); vertVc= zeros(2, 3);
idxF_new= zeros(1, length(idxF));
endP_A=  NaN*zeros(length(idxF), 3);     endP_B=  NaN*zeros(length(idxF), 3);
endV_A=  NaN*zeros(length(idxF), 3);     endV_B=  NaN*zeros(length(idxF), 3);
dirVc_A= NaN*zeros(length(idxF), 3);     dirVc_B= NaN*zeros(length(idxF), 3);

count= 0;
if ~endFlag
    for k= 1:length(idxF)
        
        % suche ob Kurven anfang in ROI ...
        if ~endFlag
            segIdx= ftrStruct.connectCell{idxF(k)};
            if segIdx(1) > 0
                dummyVc_A= ftrStruct.curveSegCell{segIdx(1)}(1:2, :);
            else
                dummyVc_A= ftrStruct.curveSegCell{-segIdx(1)}(end:(end - 1), :);
            end
        else
            dummyVc_A= reshape(ftrStruct.user.endA_Vc(idxF(k), :), [3 2])';
        end
        tmpVc= round(0.5*(dummyVc_A(1, :) + dummyVc_A(2, :)));
        vertVc_A= round(dummyVc_A(1, :) - tmpVc) + tmpVc;
        ok_A= any(((roiVert(:, 1) == vertVc_A(1) & (roiVert(:, 2) == vertVc_A(2)) & (roiVert(:, 3) == vertVc_A(3)))'));
        %    ok_A= max(sum(roiVert == (ones(size(roiVert, 1), 1)*vertVc_A), 2)) == 3;
        if (ok_A && ~bothFlag)
            count= count + 1;
            idxF_new(count)= k;
            endP_A(count, :)= dummyVc_A(1, :);
            endV_A(count, :)= vertVc_A(1, :);
            dirVc_A(count, :)= dummyVc_A(1, :) - dummyVc_A(2, :);
        elseif (ok_A == bothFlag) % falls nicht suche ob Kurven ende in ROI
            if ~endFlag
                if segIdx(end) > 0
                    dummyVc_B= ftrStruct.curveSegCell{segIdx(end)}((end - 1):end, :);
                else
                    dummyVc_B= ftrStruct.curveSegCell{-segIdx(end)}(2:-1:1, :);
                end
            else
                dummyVc_B= reshape(ftrStruct.user.endB_Vc(idxF(k), [4 5 6 1 2 3]), [3 2])';
            end
            tmpVc= round(0.5*(dummyVc_B(1, :) + dummyVc_B(2, :)));
            vertVc_B= round(dummyVc_B(2, :) - tmpVc) + tmpVc;
            ok_B= any(((roiVert(:, 1) == vertVc_B(1) & (roiVert(:, 2) == vertVc_B(2)) & (roiVert(:, 3) == vertVc_B(3)))'));
            if ok_B
                count= count + 1;
                idxF_new(count)= k;
                if bothFlag == 1
                    endP_A(count, :)= dummyVc_A(1, :);              endV_A(count, :)= vertVc_A(1, :);
                    dirVc_A(count, :)= dummyVc_A(1, :) - dummyVc_A(2, :);
                    endP_B(count, :)= dummyVc_B(2, :);              endV_B(count, :)= vertVc_B(1, :);
                    dirVc_B(count, :)= dummyVc_B(2, :) - dummyVc_B(1, :);
                else
                    endP_A(count, :)= dummyVc_B(2, :);              endV_A(count, :)= vertVc_B(1, :);
                    dirVc_A(count, :)= dummyVc_B(2, :) - dummyVc_B(1, :);
                end                
            end
        end
        if verbose && (mod(k, 100) == 0)
            verbStr= sprintf('Comparing: %3.1f%% curves finished ', round(1000*k/length(idxF))/10);
            if verbose
%                 set(txtHandle, 'String', verbStr);
%                 drawnow;
            else
                if (mod(k, 4000) == 0)
                    disp(verbStr);
                end
            end
        end
    end
else % if the endpoints are saved separately ....
    roiNo= size(roiVert, 1);
    packSizeMax= ceil(3000/roiNo);  % Seems to be near the optimum for speed

    if length(idxF) < packSizeMax
        packSize= length(idxF);
    elseif 2 > packSizeMax
        packSize= 2;
    else
        packSize= packSizeMax;
    end
    packNo= ceil(length(idxF)/packSize);

    if verbose
        set(txtHandle, 'String', sprintf('ftrstruct_query::endPointInArea: Prepare data (pack size= %d by %d)', packSize, packNo));
        drawnow;
    end

    roiMesh= zeros(packSize, roiNo, 3);
    for k= 1:3
        roiMesh(:, :, k)= roiVert(:, k*ones(packSize, 1))';
    end
    tmpMesh= zeros(packSize, roiNo);
    cmpA= zeros(packSize, 1);
    cmpB= zeros(packSize, 1);

    for i= 1:packNo
        idxPackF= idxF(floor(1 + (i - 1)*length(idxF)/packNo):floor(i*length(idxF)/packNo));
        
        % Calc endPoints A
        end_Last= ftrStruct.user.endA_Vc(idxPackF, [1 2 3]);
        end_BeforeLast= ftrStruct.user.endA_Vc(idxPackF, [4 5 6]);
        tmpVc= round(0.5*(end_BeforeLast + end_Last));
        endA= round(end_Last - tmpVc) + tmpVc;
        endPA= end_Last;
        dirA= end_Last - end_BeforeLast;
        
        % Calc endPoints B
        end_Last= ftrStruct.user.endB_Vc(idxPackF, [1 2 3]);
        end_BeforeLast= ftrStruct.user.endB_Vc(idxPackF, [4 5 6]);
        tmpVc= round(0.5*(end_BeforeLast + end_Last));
        endB= round(end_Last - tmpVc) + tmpVc;
        endPB= end_Last;
        dirB= end_Last - end_BeforeLast;
        
        % Starte Vergleich
        cmpA= max((endA(:, 1*ones(roiNo, 1)) == roiMesh(1:length(idxPackF), :, 1)) ...
            & (endA(:, 2*ones(roiNo, 1)) == roiMesh(1:length(idxPackF), :, 2)) ...
            & (endA(:, 3*ones(roiNo, 1)) == roiMesh(1:length(idxPackF), :, 3)), [], 2);
        
        cmpB= max((endB(:, 1*ones(roiNo, 1)) == roiMesh(1:length(idxPackF), :, 1)) ...
            & (endB(:, 2*ones(roiNo, 1)) == roiMesh(1:length(idxPackF), :, 2)) ...
            & (endB(:, 3*ones(roiNo, 1)) == roiMesh(1:length(idxPackF), :, 3)), [], 2);

%       [max(endA(:, 3)) max(endB(:, 3)) min(roiVert(:, 3)) max(roiVert(:, 3)) sum(sum(cmpA)) sum(sum(cmpB)) sum(sum(cmpA & cmpB))]
        if bothFlag
            fibIdx= find(cmpA & cmpB);
            endV_A((count + 1):(count + length(fibIdx)), :)= endA(fibIdx, :);
            endP_A((count + 1):(count + length(fibIdx)), :)= endPA(fibIdx, :);
            dirVc_A((count + 1):(count + length(fibIdx)), :)= dirA(fibIdx, :);
            endV_B((count + 1):(count + length(fibIdx)), :)= endB(fibIdx, :);
            endP_B((count + 1):(count + length(fibIdx)), :)= endPB(fibIdx, :);
            dirVc_B((count + 1):(count + length(fibIdx)), :)= dirB(fibIdx, :);
            idxF_new((count + 1):(count + length(fibIdx)))= floor((i - 1)*length(idxF)/packNo) + fibIdx;
            count= count + length(fibIdx);
        else
            fibIdx= find(cmpA); % side A
            endV_A((count + 1):(count + length(fibIdx)), :)= endA(fibIdx, :);
            endP_A((count + 1):(count + length(fibIdx)), :)= endPA(fibIdx, :);
            dirVc_A((count + 1):(count + length(fibIdx)), :)= dirA(fibIdx, :);
            idxF_new((count + 1):(count + length(fibIdx)))= floor((i - 1)*length(idxF)/packNo) + fibIdx;
            count= count + length(fibIdx);

            fibIdx= find(cmpB); % side B
            endV_A((count + 1):(count + length(fibIdx)), :)= endB(fibIdx, :);
            endP_A((count + 1):(count + length(fibIdx)), :)= endPB(fibIdx, :);
            dirVc_A((count + 1):(count + length(fibIdx)), :)= dirB(fibIdx, :);
            idxF_new((count + 1):(count + length(fibIdx)))= floor((i - 1)*length(idxF)/packNo) + fibIdx;
            count= count + length(fibIdx);
        end
        if verbose && (mod(i, floor(packNo/20)) == 0)
            verbStr= sprintf('Comparing: %3.1f%% curves finished ', round(1000*i/packNo)/10);
            if verbose
                set(txtHandle, 'String', verbStr);
                drawnow;
            else
                if (mod(k, 4000) == 0)
                    disp(verbStr);
                end
            end
        end
    end
end


if count == 0
    fiberStruct= [];
    errStr= 'ftrstruct_query (local_getInArea): No fibers found';
else
    fiberStruct= ftrstruct_init('fiber');
    fiberStruct.curveID= idxF(idxF_new(1:count));
    if bothFlag
        fiberStruct.user.endP_A=    endP_A(1:count, :);        
        fiberStruct.user.endV_A=    endV_A(1:count, :);        
        fiberStruct.user.dirVc_A=  dirVc_A(1:count, :);
        fiberStruct.user.endP_B=    endP_B(1:count, :);        
        fiberStruct.user.endV_B=    endV_B(1:count, :);        
        fiberStruct.user.dirVc_B=  dirVc_B(1:count, :);
    elseif isempty(nameRef)
        fiberStruct.user.endP_A=    endP_A(1:count, :);        
        fiberStruct.user.endV_A=    endV_A(1:count, :);        
        fiberStruct.user.dirVc_A=  dirVc_A(1:count, :);
    else    
        fiberOrg= [];
        if ischar(nameRef)
            if not(isempty(private_name2id(ftrStruct, nameRef)));
                fiberOrg= ftrStruct.fiber{private_name2id(ftrStruct, nameRef)};
            end
        end;
        if ~isempty(fiberOrg) && isfield(fiberOrg.user, 'endP_A') && isfield(fiberOrg.user, 'endV_A')
            fiberStruct.user.endP_A= fiberOrg.user.endP_A(idxF_new(1:count), :);
            fiberStruct.user.endV_A= fiberOrg.user.endV_A(idxF_new(1:count), :);
            fiberStruct.user.dirVc_A=  fiberOrg.user.dirVc_A(idxF_new(1:count), :);
            fiberStruct.user.endP_B= endP_A(1:count, :);        
            fiberStruct.user.endV_B= endV_A(1:count, :);        
            fiberStruct.user.dirVc_B=dirVc_A(1:count, :);
        else
            fiberStruct.user.endP_A= endP_A(1:count, :);        
            fiberStruct.user.endV_A= endV_A(1:count, :);        
            fiberStruct.user.dirVc_A=dirVc_A(1:count, :);
        end
    end
end

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
%
%  START:
%       [id, errStr]= private_name2id(ftrStruct, nameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id, errStr]= private_name2id(ftrStruct, nameStr)
errStr= '';

if isnumeric(nameStr)
    id= nameStr;
    if (id < 1) || (id > length(ftrStruct.fiber))
        id= [];
        errStr= 'ftrstruct_query (private_name2id): invalide fiber id';
    end
    return
end

id= [];
for i= 1:length(ftrStruct.fiber)
    if strcmp(ftrStruct.fiber{i}.name, nameStr)
        id= i;
        return
    end
end

if isempty(id)
    errStr= strcat('ftrstruct_query (private_name2id): ''', nameStr, ''' is an invalide fiber name');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: private_name2id
%

%
%
%  START:
%       [id, errStr]= private_getFiberId(ftrStruct, nameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idxF, errStr, fstr]= private_getFiberId(ftrStruct, nameStr)
errStr= ''; idxF= []; fstr = [];
if isempty(nameStr) || strcmp(nameStr, 'all')
    idxF= 1:length(ftrStruct.connectCell);
    return;
end
id= private_name2id(ftrStruct, nameStr);
if isempty(id)
    errStr= sprintf('Error(ftrstruct_query::private_getFiberId): Fiber ''%s'' is unknown', nameStr);
    return
end
idxF= ftrStruct.fiber{id}.curveID;
fstr= ftrStruct.fiber{id};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: private_name2id
%

%
%
%  START:
%       [res, errStr]= private_getVertex(ftrStruct, ratio)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_getVertex(conCell, vertCell, id, type, trans2)

errStr= ''; res= [];

if ~exist('trans2') || isempty(trans2)
    trans2= [];
end

lenCrv= zeros(length(id), 1);
resCell= {};
for j= 1:length(id)
    idx= conCell{id(j)};    len= zeros(length(idx), 1);
    for i= 1:length(idx)
        len(i)= size(vertCell{abs(idx(i))}, 1);
    end
    resCell{j, 1}= zeros(sum(len), 3);
    pos= 1;
    for i= 1:length(idx)
        if idx(i) > 0
            resCell{j, 1}(pos:(pos + len(i) - 1), :)= vertCell{idx(i)};
        else
            vertTemp= vertCell{-idx(i)};
            resCell{j, 1}(pos:(pos + len(i) - 1), :)= vertTemp(end:-1:1, :);
        end
        pos= pos + len(i);
    end
    if ~isempty(trans2)
        resCell{j, 1}= resCell{j, 1}*(trans2(1:3, 1:3)');
%        resCell{j, 1}= trans2(1:3, 1:3)*reshape(resCell{j, 1}, [2 1]);
        resCell{j, 1}(:, 1)= resCell{j, 1}(:, 1) + trans2(1, 4);
        resCell{j, 1}(:, 2)= resCell{j, 1}(:, 2) + trans2(2, 4);
        resCell{j, 1}(:, 3)= resCell{j, 1}(:, 3) + trans2(3, 4);
    end
    lenCrv(j)= sum(len);
end

if strcmp(type, 'array')
    res= zeros(sum(lenCrv) + length(id) - 1, 3);
    pos= 1;
    for i= 1:length(id)
        res(pos:(pos + lenCrv(i)), :)= [resCell{i}; [NaN NaN NaN]];
        pos= pos + lenCrv(i) + 1;
    end
elseif strcmp(type, 'cell')
    res= resCell;
else
    errStr= 'ftrstruct_query (private_getVertex): inavlid type';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: private_getVertex
%

%
%
%  START:
%       [res, errStr]= private_getHMatrix(ftrStruct, flag)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_getHMatrix(ftrStruct, flag)

res= []; errStr= '';

if isempty(flag) || (flag == 'm') || (flag == 'M')
    return
elseif isequal(size(flag), [4 4])
    res= flag;
elseif isequal(size(flag), [3 3])
    res= diag(ones(4, 1));
    res(1:3, 1:3)= flag;
elseif isequal(size(flag), [3 4])
    res= diag(ones(4, 1));
    res(1:3, 1:4)= flag;
elseif (flag == 'v') || (flag == 'V')
    res= diag([ftrStruct.vox(1) ftrStruct.vox(2) ftrStruct.vox(3) + ftrStruct.vox(4) 1]);
elseif (flag == 'w') || (flag == 'W')
    [verStr, errStr]= ftrstruct_query(ftrStruct, 'getVer');
    if ~strcmp('V1.0', verStr)
        if ~isequal(size(ftrStruct.hMatrix), [4 4]) || ~isnumeric(ftrStruct.hMatrix)
            errStr= strcat(mfilename, '::private_getHMatrix (error): hMatrix entry is not a 4x4 matrix. Using identity');
        else
            res= ftrStruct.hMatrix;
        end
    else
        errStr= strcat(mfilename, '::private_getHMatrix (error): Version of ftrStruct (V1.0) does not support transformation');
        return
    end
else
    errStr= strcat(mfilename, '::private_getHMatrix (error): Unrecognized flag for transformation');
    return
end