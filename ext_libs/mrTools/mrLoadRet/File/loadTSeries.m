function tseries = loadTSeries(view,scan,slice,frame,x,y,precision, groupNum)
%
% tSeries = loadTSeries(view,[scan],[slice],[frame],[x],[y],[precision], [groupNum])
%
% Loads the tSeries corresponding to the specified scan and slice. The
% tseries files must be in <homedir>/<view>/<group>/TSeries, and must be
% nifti (or analyze) format (4D file containing all of the slices and
% temporal frames from a given scan).
%
% scan: number specfying which scan to load.
%
% slice: either a number specifying which slice to load or 'all' which
% loads all of the slices. Default: 'all'
%
% frame: optional number specifying which temporal frame to load. Ignored
% if slice is specified. Default: [] (load all frames).
%
% precision: precision of data returned, default is 'double' can be 'single'
%
% tseries: If a single slice is specified, tseries is returned as aa 2D
% matrix (nFrames x nVoxels). If slice is 'all' then the tseries is
% returned as a 4D array [x y z t]. If frame is specified, then the tseries
% is a single volume returned as a 3D array [x y t]
%
% djh 1/9/98
% djh 2/20/2001 Removed interpolation, dumped dat files
% JL 8/6/03 Updated to allow Analyze format.
% djh 5/2005 Updated to mrLoadRet 4.0
% jlg 3/2007 Streamlined code, also selection of a set of
%            slices/frames by adding argument slice=[min max], or
%            frame=[min max]. No longer reshapes tseries for single slice
% JB 07/2011 Added groupNum option

if ieNotDefined('scan')
  scan = viewGet(view,'curScan');
end
if ieNotDefined('slice')
  slice = 'all';
end
if ieNotDefined('frame')
  frame = [];
end
if ieNotDefined('x')
  x = [];
end
if ieNotDefined('y')
  y = [];
end
if ieNotDefined('precision'), precision = mrGetPref('defaultPrecision');end
if ieNotDefined('groupNum'), groupNum = viewGet(view,'curGroup');end

% Get the pathStr to the tseries file
pathStr = viewGet(view,'tseriesPathStr',scan,groupNum);
% Error if file not found
if ~exist(pathStr,'file')
  mrErrorDlg(['File ',pathStr,' not found']);
end

% change slices to an array
if strcmp(slice,'all'),slice = [];,end
% and validate it
if ~isempty(slice) & ((1 > min(slice)) | (max(slice) > viewGet(view,'nslices',scan,groupNum)))
  mrErrorDlg(['Invalid slice number: ',num2str(slice)]);
  keyboard
end

% validate frames
if ~isempty(frame) & ((1 > frame) | (frame > viewGet(view,'nframes',scan,groupNum)))
  mrErrorDlg(['Invalid frame number: ',num2str(frame)]);
end

% Load it
[tseries,hdr] = mlrImageReadNifti(pathStr,{x,y,slice,frame},precision);
dims = size(tseries);

% check frame count
if isempty(frame)
  nFramesExpected = viewGet(view,'totalFrames',scan,groupNum);
elseif length(frame) == 1
  nFramesExpected = 1;
elseif length(frame) == 2
  nFramesExpected= frame(2)-frame(1)+1;
end
% nFrames = dims(4);
% if (nFrames ~= nFramesExpected)
%   mrWarnDlg('loadTSeries: number of frames in tseries file does not match expected.');
% end
