function xyzmm = ea_getplanes(varargin)
% Get current slice indices or millimeter coordinates from 3D scence.
%
% By default, the output is millimeter coordinates

outputtype = 'mm';

if nargin < 1
    resultfig = gcf;
elseif ishandle(varargin{1})
    resultfig = varargin{1};
    if nargin == 2
        outputtype = varargin{2};
    end
elseif ischar(varargin{1})
    resultfig = gcf;
    outputtype = varargin{1};
end

xslice = getappdata(resultfig, 'xsliceplot');
yslice = getappdata(resultfig, 'ysliceplot');
zslice = getappdata(resultfig, 'zsliceplot');

xyzmm = [xslice.UserData.sliceidx, yslice.UserData.sliceidx, zslice.UserData.sliceidx];

if strcmp(outputtype, 'mm')
    xyzmm = ea_vox2mm(xyzmm, xslice.UserData.I2X);
end
