function [res, errStr]= morph_data_ext(varargin)
%function [res, errStr]= morph_data_ext(mask, commandStr, par1, par2, par3)
%
% Funktion to apply any binary morphology kernel an binary images or 3D datasets.
%
% Input: 
%       -mask: mrStruct and schould contain a 2D or 3D data. All values ~=
%              zeros are set to one in mask
%       -kernelAy: Must have the same dimension as mask. The size of
%                  dimensions a free. 
%       -kernelPos: Describes the correspondence between mask and kernelAy
%       -OPStr: {'<='; '>='; '<'; '>'; '=='; '~='} it describes the
%               condition if a voxel will be set or not
%       -noVoxelInKernel: the number on the right side of OPStr
%       -forceKPos: {'yes'; 'no'} if set 'yes', the voxel X in res can only
%                   be set if it is also set in mask
%
% 10/06
%
% UNIX


res= []; errStr= [];

if nargin < 2
    errStr= strcat(mfilename, '(error): funktion needs at least two parameter');
    return
end

if ~mrstruct_istype(varargin{1}) && (length(mrstruct_query(varargin{1}, 'sizeAy')) == 3)
    errStr= strcat(mfilename, '(error): first parameter have to be a volume mrstruct');
    return
end
mrStruct= varargin{1};

if ~ischar(varargin{2})
    errStr= strcat(mfilename, '(error): second parameter have to be a string');
    return;
end
commandStr= varargin{2};

for i= 1:7
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end



if strcmp(commandStr, 'growByCircle')
    rad= param{1};
    while isempty(rad)
        prompt={'Please enter the radius of the growing sphere [mm]:'};
        def={'10'};
        dlgTitle='Growing Sphere';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        if isempty(answer) 
            errStr= strcat(mfilename, '(error): Aborted by user');
            return;            
        end
        rad= str2num(answer{1});
    end
    if ishandle(param{2})
        verHd= param{2};
    else 
        verHd= [];
    end
    [res, errStr]= local_growByCircle(mrStruct, rad, verHd);
    return;
else
    errStr= sprintf('%s(error): command ''%s'' is not implemented', mfilename, commandStr);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [res, errStr]= local_growByCircle(mrStruct, rad, verHd)

res= []; errStr= '';

sizeAy= mrstruct_query(mrStruct, 'sizeAy');

messStr= sprintf('%s: Prepare data', mfilename);
if isempty(verHd)
    disp(messStr)
else
    set(verHd, 'String', messStr);
    drawnow;
end

% berechner randpixel
kernelAy= zeros(3, 3, 3);
kernelAy(:, 2, 2)= 1; kernelAy(2, :, 2)= 1; kernelAy(2, 2, :)= 1;
tmp= morph_data(mrStruct, kernelAy, [2 2 2], '<', 7, 'yes');
idx= find(tmp.dataAy ~= 0);

seedVc= reshape_index(idx, sizeAy);
sphereVc= privat_genSphere(mrStruct.vox, rad);
sphLen= size(sphereVc, 1);

res= mrStruct;

for i= 1:size(seedVc, 1)
    markVc= sphereVc + ones(sphLen, 1)*seedVc(i, :);    
    for j= 1:3
        idx= find((markVc(:, j) < 1) | (markVc(:, j) > sizeAy(j)));
        if ~isempty(idx)
            markVc(idx, :)= [];
        end
    end

    messStr= sprintf('%s: Processing ... %d/%d', mfilename, i, size(seedVc, 1));
    if isempty(verHd)
        disp(messStr)
    else
        set(verHd, 'String', messStr);
        drawnow;
    end
    
    markIdx= reshape_index_back(markVc, sizeAy);
    res.dataAy(markIdx)= 1;
end

   messStr= sprintf('%s: Finished', mfilename);
    if isempty(verHd)
        disp(messStr)
    else
        set(verHd, 'String', messStr);
        drawnow;
    end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vertVc= privat_genSphere(voxSize, rad)

if isempty(voxSize)
    voxSize= ones(1, 3);
    warning('morph_data_ext:privat_genSphere(warning) No voxel size is given. Assume [1 1 1 0]');
end
pNo= ceil(rad./voxSize(1:3));
[X, Y, Z]= meshgrid((0:pNo(1))*voxSize(1), (0:pNo(2))*voxSize(2), (0:pNo(3))*voxSize(3));
idx= find(sqrt(X.^2 + Y.^2 + Z.^2) <= rad);
tmpVc= reshape_index(idx, size(X)) - 1;


vertVc= [...
        [ tmpVc(:, 1),  tmpVc(:, 2),  tmpVc(:, 3)]; ...
        [ tmpVc(:, 1),  tmpVc(:, 2), -tmpVc(:, 3)]; ...
        [ tmpVc(:, 1), -tmpVc(:, 2),  tmpVc(:, 3)]; ...
        [ tmpVc(:, 1), -tmpVc(:, 2), -tmpVc(:, 3)]; ...
        [-tmpVc(:, 1),  tmpVc(:, 2),  tmpVc(:, 3)]; ...
        [-tmpVc(:, 1),  tmpVc(:, 2), -tmpVc(:, 3)]; ...
        [-tmpVc(:, 1), -tmpVc(:, 2),  tmpVc(:, 3)]; ...
        [-tmpVc(:, 1), -tmpVc(:, 2), -tmpVc(:, 3)]];
