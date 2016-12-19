%function [res, errStr]= mrtsruct_transfortm(mrStruct, commandStr[, op1[, op2[... [, opN]]]])
%
%   command:    {'x-axis' | 'y-axis' | 'z-axis | 'translate'}
%   
%
% Bjoern W. Kreher
% 04/04
%
% UNIX

function [res, errStr]= mrstruct_transform(varargin)

res= []; errStr= '';
maxArgs= 7;

%% Parameter Kontrolle
mrStruct= varargin{1};
commStr= varargin{2};
if ~mrstruct_istype(mrStruct)
    errStr= 'mrstruct_transform: First parametret is not a valid mrstruct';
    return;
end
if (mrstruct_query(mrStruct, 'dimensions') < 2) || (mrstruct_query(mrStruct, 'dimensions') > 3)
    errStr= 'mrstruct_transform: mrStruct has to contain at least two dimensions';
    return
end
if ~ischar(commStr)
    errStr= 'mrstruct_transform: Second parametre have to be a string';
    return;
end

% Uebergabe dynamischer parameter
argsCell= cell(maxArgs, 1);
for i= 1:maxArgs
    if i <= (nargin - 2)
        argsCell{i, 1}= varargin{i + 2};
    else
        argsCell{i, 1}= [];
    end
end

if strcmp(commStr, 'x-axis')
    angle= argsCell{1, 1};  centerVc= argsCell{2, 1};
    [res, errStr]= local_rotate_X(mrStruct, angle, centerVc);
elseif strcmp(commStr, 'y-axis')
    angle= argsCell{1, 1};  centerVc= argsCell{2, 1};
    [res, errStr]= local_rotate_Y(mrStruct, angle, centerVc);
elseif strcmp(commStr, 'z-axis')
    angle= argsCell{1, 1};  centerVc= argsCell{2, 1};
    [res, errStr]= local_rotate_Z(mrStruct, angle, centerVc);
    
end



%%%
%%% START: function [res, errStr]= local_rotate_X(mrStruct, angle, centerVc)
function [res, errStr]= local_rotate_X(mrStruct, angle, centerVc)

res= []; errStr= '';

if isempty(centerVc)
    [centerVc, errStr]= private_calcMassCenter(mrStruct);
    if isempty(centerVc)
        return
    end
end

sizeAy= mrstruct_query(mrStruct, 'sizeAy');

for i= 1:sizeAy(1)
    mrStruct.dataAy(i, :, :)= private_rotateSlice(squeeze(mrStruct.dataAy(i, :, :)), angle, centerVc([2 3]), mrStruct.vox([2 3]));
    if mod(i, 10) == 0
        fprintf('%d/%d Schichten rotiert', i, sizeAy(1));
    end
end

if isfield(mrStruct.user, 'transform')
    invert= mrStruct.user.transform.invert;
    hist= mrStruct.user.transform.hist;
    rotMx= mrStruct.user.transform.rotMx;
else
    invert= 'mrStruct';
    hist= 'mrStruct';
    rotMx= [];
end

idx= strfind(invert, 'mrStruct');
mrStruct.user.transform.invert= sprintf('%smrstruct_transform(mrStruct, ''x-axis'', %d, [%s])%s', invert(1:(idx-1)), -angle, vector2str(centerVc), invert((idx + 8):end));
mrStruct.user.transform.hist= sprintf('mrstruct_transform(%s, ''x-axis'', %d, [%s])', invert, angle, vector2str(centerVc));
mrStruct.user.transform.rotMx= hm_trans(centerVc, hm_rot(angle, 'x-axis', hm_trans(-centerVc, rotMx)));
res= mrStruct;

%%% END: local_rotate_X
%%%

%%%
%%% START: function [res, errStr]= local_rotate_Y(mrStruct, angle, centerVc)
function [res, errStr]= local_rotate_Y(mrStruct, angle, centerVc)

res= []; errStr= '';

if isempty(centerVc)
    [centerVc, errStr]= private_calcMassCenter(mrStruct);
    if isempty(centerVc)
        return
    end
end

sizeAy= mrstruct_query(mrStruct, 'sizeAy');

for i= 1:sizeAy(2)
    mrStruct.dataAy(:, i, :)= private_rotateSlice(squeeze(mrStruct.dataAy(:, i, :)), angle, centerVc([1 3]), mrStruct.vox([1 3]));
    if mod(i, 10) == 0
        fprintf('%d/%d Schichten rotiert', i, sizeAy(2));
    end
end

if isfield(mrStruct.user, 'transform')
    invert= mrStruct.user.transform.invert;
    hist= mrStruct.user.transform.hist;
    rotMx= mrStruct.user.transform.rotMx;
else
    invert= 'mrStruct';
    hist= 'mrStruct';
    rotMx= [];
end

idx= strfind(invert, 'mrStruct');
mrStruct.user.transform.invert= sprintf('%smrstruct_transform(mrStruct, ''y-axis'', %d, [%s])%s', invert(1:(idx-1)), -angle, vector2str(centerVc), invert((idx + 8):end));
mrStruct.user.transform.hist= sprintf('mrstruct_transform(%s, ''y-axis'', %d, [%s])', invert, angle, vector2str(centerVc));
mrStruct.user.transform.rotMx= hm_trans(centerVc, hm_rot(angle, 'y-axis', hm_trans(-centerVc, rotMx)));
res= mrStruct;

%%% END: local_rotate_Y
%%%

%%%
%%% START: function [res, errStr]= local_rotate_Z(mrStruct, angle, centerVc)
function [res, errStr]= local_rotate_Z(mrStruct, angle, centerVc)

res= []; errStr= '';

if isempty(centerVc)
    [centerVc, errStr]= private_calcMassCenter(mrStruct);
    if isempty(centerVc)
        return
    end
end

sizeAy= mrstruct_query(mrStruct, 'sizeAy');

for i= 1:sizeAy(3)
    mrStruct.dataAy(:, :, i)= private_rotateSlice(mrStruct.dataAy(:, :, i), angle, centerVc([1 2]), mrStruct.vox([1 2]));
    if mod(i, 10) == 0
        fprintf('%d/%d Schichten rotiert', i, sizeAy(3));
    end
end

if isfield(mrStruct.user, 'transform')
    invert= mrStruct.user.transform.invert;
    hist= mrStruct.user.transform.hist;
    rotMx= mrStruct.user.transform.rotMx;
else
    invert= 'mrStruct';
    hist= 'mrStruct';
    rotMx= [];
end

idx= strfind(invert, 'mrStruct');
mrStruct.user.transform.invert= sprintf('%smrstruct_transform(mrStruct, ''z-axis'', %d, [%s])%s', invert(1:(idx-1)), -angle, vector2str(centerVc), invert((idx + 8):end));
mrStruct.user.transform.hist= sprintf('mrstruct_transform(%s, ''z-axis'', %d, [%s])', invert, angle, vector2str(centerVc));
mrStruct.user.transform.rotMx= hm_trans(centerVc, hm_rot(angle, 'z-axis', hm_trans(-centerVc, rotMx)));
res= mrStruct;

%%% END: local_rotate_Z
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%% START: function [res, errStr]= private_rotateSlice(sliceAy, angle, centerVc)
function [res, errStr]= private_rotateSlice(sliceMy, angle, centerVc, vox)

res= []; errStr= '';

sizeAy= size(sliceMy);

% Konten an alle Voxel
vertVc= reshape_index(1:prod(sizeAy), sizeAy);
% Translation in das Rotations Zentrum
vertVc= vertVc - (ones(prod(sizeAy), 1)*centerVc);
% Scalierung auf voxelgroesse
vertVc= vertVc.*(ones(prod(sizeAy), 1)*vox);
% Nagative Rotatio der Knoten
vertVc= [vertVc(:, 1)*cos(angle) + vertVc(:, 2)*sin(angle) -vertVc(:, 1)*sin(angle) + vertVc(:, 2)*cos(angle)];
% Scalierung zurck in voxelabstnden
vertVc= vertVc.*(ones(prod(sizeAy), 1)*(1./vox));
% Translation in Ursprung
vertVc= vertVc + (ones(prod(sizeAy), 1)*centerVc);
% Eleminierung von Knoten ausserahlb im Bereich

vertA_Vc= private_eliminate_outside([floor(vertVc(:, 1)) floor(vertVc(:, 2))], sizeAy);
vertB_Vc= private_eliminate_outside([ceil(vertVc(:, 1))  floor(vertVc(:, 2))], sizeAy);
vertC_Vc= private_eliminate_outside([floor(vertVc(:, 1)) ceil(vertVc(:, 2))], sizeAy);
vertD_Vc= private_eliminate_outside([ceil(vertVc(:, 1))  ceil(vertVc(:, 2))], sizeAy);

xLin= vertVc(:, 1) - floor(vertVc(:, 1)); 
yLin= vertVc(:, 2) - floor(vertVc(:, 2)); 

idxA= reshape_index_back(round(vertA_Vc), sizeAy);
idxB= reshape_index_back(round(vertB_Vc), sizeAy);
idxC= reshape_index_back(round(vertC_Vc), sizeAy);
idxD= reshape_index_back(round(vertD_Vc), sizeAy);

res= reshape((sliceMy(idxA).*(1 - xLin) + sliceMy(idxB).*xLin).*(1 - yLin) + (sliceMy(idxC).*(1 - xLin) + sliceMy(idxD).*xLin).*yLin, sizeAy);

%%% END: private_rotateSlice
%%%


%%%
%%% START: function [res, errStr]= private_calcMassCenter(mrStruct)
function [res, errStr]= private_calcMassCenter(mrStruct)

res= []; errStr= '';
if isempty(mrStruct.dataAy)
    errStr= 'mrstruct_transfor: mrStruct is empty';
    return
end

% maske zur schwerpunktsberechnung entspricht allen voxeln die groesser als mittelwert/2 sind
sizeAy= mrstruct_query(mrStruct, 'sizeAy');
meanVal= mean(reshape(mrStruct.dataAy, [1, prod(sizeAy)]));
res= mean(reshape_index(find(mrStruct.dataAy > meanVal/2), sizeAy));

%%% END private_calcMassCenter
%%%

%%%
%%% START: function res= eliminate_outside(vectsVc, sizeAy)
function vectsVc= private_eliminate_outside(vectsVc, sizeAy, rep)

if ~exist('rep') || isempty(rep)
    rep= ones(1, length(sizeAy));
end

for i= 1:length(sizeAy)
    idx= vectsVc(:, i) < 1;
%    vectsVc(idx, :)= ones(length(idx), 1)*rep;
    vectsVc(idx, :)= 1;
    idx= vectsVc(:, i) > sizeAy(i);
%    vectsVc(idx, :)= ones(length(idx), 1)*rep;
    vectsVc(idx, :)= sizeAy(i);
end
%%% END: eliminate_outside
%%%