% mrSliceExport - output series of images across all slices
%
%      usage: [  ] = mrSliceExportTest( <thisView>, <imrange>, <sliceList>, <fname>, <nRows>, <nColumns>, <scanList>) )
%         by: denis schluppeck, modified by julien besle for integration in mrLoadRet
%       date: 2007-11-22
%        $Id: mrSliceExportCor.m 583 2010-10-15 09:54:53Z lpzjb $:
%     inputs:       thisView: mrLoadRet View structure
%             imrange: relative dimensions of the display to export between 0 and 1 ([minX maxX minY maxY] form left to right and top to bottom)
%           sliceList:
%               fname: base name of the tif file to export to
%               nRows: number of rows in the image
%            scanList:
%    outputs: 
%
%    purpose: exports currently displayed base+overlay as a tif image, (goes across required scans and slices)
%
%        e.g:
%
function [  ]=mrSliceExport( thisView, imrange, sliceList, fname, nRows, scanList)

mrGlobals

if ieNotDefined('imrange')
  imrange = [0 1 0 1];
end
% clip between 0 and 1
imrange = min(imrange,1);
imrange = max(imrange,0);
  
if ieNotDefined('thisView')
  thisView = viewGet([],'view',viewGet([],'viewnums'));
  disp('... using current view')
  %disp(thisView)
end

if ieNotDefined('fname')
   fname = 'montage00.tif';
end

baseType=viewGet(thisView,'basetype');
switch(baseType)
  case 0
    if ieNotDefined('sliceList')
      sliceList = viewGet(thisView,'nslices'):-1:1;% 1:nSlices;
    end
  case{1,2}
    if ieNotDefined('sliceList')
     %if sliceList is empty and we have a flat or a surface base, we just want to export the image in the current settings of depths
      baseType= -1;
      sliceList=1;
    end
end

if ieNotDefined('scanList')
   scanList = viewGet(thisView,'curscan');
else
   nScans = viewGet(thisView,'nscans');
   if any(scanList>nScans)
      mrWarnDlg(['There are only ' num2str(nScans) ' scans']);
      return;
   end
end

nSlices = length(sliceList);
if ieNotDefined('nRows')
  nRows = 2;
end
nColumns= ceil(length(scanList)*nSlices/nRows);




if isempty(viewGet(thisView,'curanalysis'))
  % load in current analaysis
  thisView = loadAnalysis(thisView, 'corAnal');
end

if isempty(viewGet(thisView,'curbase'))
  % load anatomy
  thisView = loadAnat(thisView);
end

% don't want to hardcode any parameters here, but you could roll your own
% e.g.
% thisView = viewSet(thisView, 'curoverlay', 1); % co is 1
% thisView = viewSet(thisView,'overlaymin', 0.35);
% thisView = viewSet(thisView,'overlaymax', 1.0); 

viewNum = viewGet(thisView,'viewnum');
curGroup = viewGet(thisView,'curgroup');
curOverlay = viewGet(thisView,'currentoverlay');

% for alpha transparency, have to go via mlrGuiSet -- is this bad?
% mlrGuiSet(thisView,'alpha', 0.75);

% loop over slices and get images

basedims = viewGet(thisView,'basedims');
newx = ceil(imrange(1:2)*(basedims(1)-1))+1; % indeces are from 1 to basedims
newy =  ceil(imrange(3:4)*(basedims(2)-1))+1;

fprintf(1,'Image coordinates: X: %d -> %d - Y: %d -> %d \n',newx(1),newx(2),newy(1),newy(2));

[pathname, fname, extension] = fileparts(fname);
fname = [fname extension];

% check that the directory exists
if ~exist(pathname, 'dir')
  mkdir(pathname)
  disp(['(mrSliceExport) Created directory ' pathname]);
else
  % rm all the im*.tif files in Etc
  unix(['rm ' pathname '/im*.tif']);
end

curscan = viewGet(thisView,'curscan');
switch baseType
  case 0
    curSlice = viewGet(thisView,'curslice');
  otherwise
    minDepth = viewGet(thisView,'corticalmindepth');
    maxDepth = viewGet(thisView,'corticalmaxdepth');
end
      
% loop through scans and slices.
for iScan = scanList
  thisView=viewSet(thisView, 'curscan', iScan);
  for iSlice = sliceList
    switch baseType
      case 0
        %set the slice
        thisView=viewSet(thisView, 'curslice', iSlice);
        disp(sprintf('rendering scan %d slice %d .',iScan,iSlice));
        imageFilename = sprintf('%s/im_%02d_%02d.tif',pathname,iScan, iSlice);
      case {1,2}
        %set the cortical depth
        mlrGuiSet(viewNum,'corticalMinDepth',iSlice);
        mlrGuiSet(viewNum,'corticalMaxDepth',iSlice);
        disp(sprintf('rendering scan %d depth %f .',iScan,iSlice));
        imageFilename = sprintf('%s/im_%02d_%f.tif',pathname,iScan, iSlice);
      case -1
        %do nothing
        disp(sprintf('rendering scan %d .',iScan));
        imageFilename = sprintf('%s/im_%02d_%f_%f.tif',pathname,iScan,minDepth, maxDepth);
    end


    img = refreshMLRDisplay(viewNum);
    reducedim = img(newy(1):newy(2),newx(1):newx(2),:);

    imwrite(double(reducedim), imageFilename, 'tif');
  end
end

%re-instal original scan and slice/depth
thisView=viewSet(thisView, 'curscan', curscan);
switch baseType
  case 0
    viewSet(thisView,'curslice',curSlice);
  otherwise
    mlrGuiSet(viewNum,'corticalMinDepth',minDepth);
    mlrGuiSet(viewNum,'corticalMaxDepth',maxDepth);
end
refreshMLRDisplay(viewNum);

% and create a montage for a final tif image (this is taken from mrCreateMontage by D. Schluppeck)
dirOutput = dir([pathname '/im*.tif']);
filenames = {dirOutput.name};
for iFile = 1:length(filenames)
  filenames{iFile} = [pathname '/' filenames{iFile}];
end

hFigure = selectGraphWin(0,'replace');
set(hFigure,'name','Exported Montage');
h_ = montage(filenames, 'size', [nRows nColumns]);
imwrite(get(h_,'cdata'), [pathname '/' fname], 'tiff');
disp(sprintf('wrote: %s', [pathname '/' fname]));

return ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range = findRange(data)

ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
    ampMin = min([ampMin min(data{scan}(:))]);
    ampMax = max([ampMax max(data{scan}(:))]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end




% mrCreateMontage - create montage from tif images in etc folder 
%
%      usage: [  ] = mrCreateMontage(  )
%         by: denis schluppeck
%       date: 2008-09-03
%        $Id$:
%     inputs: 
%    outputs: 
%
%    purpose: get tif files form Etc folder and make a montage
%
%        e.g:
%
function [  ]=mrCreateMontage( fname, nRows, nColumns  )














