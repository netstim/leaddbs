% interface program for batch program
% Operations with streamline tracts (Mori Tracts).
% Author: Susanne Schnell
% Linux 02.03.2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = dti_roiop_ui(cmd, P)

switch cmd
    case 'Invert'
        autofname = '_invers.mat';
        %load ROI      
        ROI = maskstruct_read(P.roiname{1});
        ROInames = maskstruct_query(ROI,'maskNames');
        % INVERT
        [ROI,errStr] = maskstruct_modify(ROI, 'inv',ROInames{P.mask1},P.maskname);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'ANDop'
        autofname = '_ANDop.mat';
        %load ROI
        ROI = maskstruct_read(P.roiname{1});
        ROInames = maskstruct_query(ROI,'maskNames');
        % AND
        [ROI,errStr] = maskstruct_modify(ROI, 'AND',ROInames{P.mask1},ROInames{P.mask2,:},P.maskname);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'ORop'
        autofname = '_ORop.mat';
        ROI = maskstruct_read(P.roiname{1});
        ROInames = maskstruct_query(ROI,'maskNames');
        % OR
        [ROI,errStr] = maskstruct_modify(ROI, 'OR',ROInames{P.mask},ROInames{P.mask2},P.maskname);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'XORop'
        autofname = '_XORop.mat';
        ROI = maskstruct_read(P.roiname{1});
        ROInames = maskstruct_query(ROI,'maskNames');
        % XOR
        [ROI,errStr] = maskstruct_modify(ROI, 'OR',ROInames{P.mask1},ROInames{P.mask2},P.maskname);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'Erosion'
        autofname = '_erosed.mat';
        % load data
        ROI = maskstruct_read(P.roiname{1});
        ROInames = maskstruct_query(ROI,'maskNames');
        % EROSION
        [ROI,errStr] = maskstruct_modify(ROI, 'erosion',ROInames{P.mask1},P.maskname,'XY',P.voxels);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'Dilation'
        autofname = '_dilated.mat';  
        % load data
        ROI = maskstruct_read(P.roiname{1});
        ROInames = maskstruct_query(ROI,'maskNames');
        % DILATION
        [ROI,errStr] = maskstruct_modify(ROI, 'dilatation',ROInames{P.mask1},P.maskname,'XY',P.voxels);
        if ~isempty(errStr)
            error(errStr)
        end

    case 'GrowSphere'
        autofname = '_sphere.mat';  
        % load data
        ROI = maskstruct_read(P.roiname{1});
        % GROW SPHERE
        [ROI,errStr] = maskstruct_modify(ROI, 'growSpheric',P.coords,P.sizeSphere{1},P.maskname);
        if ~isempty(errStr)
            error(errStr)
        end
end

% newfilename
if strcmp(P.newfilename,'.mat') || isempty(P.newfilename)
    [path,name] = fileparts(P.roiname{1});
    filename = fullfile(path, [name, autofname]);
else
    [path,name,ext] = fileparts(P.roiname{1});
    [newpath,name,ext] = fileparts(P.newfilename);
    if isempty(newpath)
        filename = fullfile(path,[name ext]);
    else
        filename = fullfile(newpath,[name ext]);
    end
end
    
%save result
maskstruct_write(ROI,filename);
out.files{1} = filename;

