% mlrImageHeaderDisp.m
%
%      usage: mlrImageDispHeader(h,<verbose=1>)
%         by: justin gardner
%       date: 09/04/11
%    purpose: displays the image header (h can be a header)
%             or can be called with a filename/view/etc. just
%             like mlrImageLoad
%
function mlrImageHeaderDisp(varargin)

% check arguments
if nargin == 0
  help mlrImageDispHeader
  return
end

% load arguments
[imageArgs otherArgs] = mlrImageParseArgs(varargin);
verbose = [];
getArgs(otherArgs,{'verbose=1'});

if length(imageArgs) < 1
  disp(sprintf('(mlrImageHeaderDisp) No files to display'));
  return
end
  
% if it is not a header already, then load it
if mlrImageIsHeader(imageArgs{1})
  h = imageArgs{1};
else
  h = mlrImageHeaderLoad(imageArgs{1});
  if isempty(h),return,end
end

% display the main part of the header
dispHeader(h.filename);
if ~isempty(h.ext)
  disp(sprintf('type: %s (%s)',h.type,h.ext));
else
  disp(sprintf('type: %s',h.type));
end
disp(sprintf('dim: [%s]',mlrnum2str(h.dim(:)','sigfigs=0')));
disp(sprintf('pixdim: [%s]',mlrnum2str(h.pixdim(:)')));
if ~isempty(h.qform)
  disp(sprintf('qform:'));
  disp(sprintf('%s',mlrnum2str(h.qform,'compact=0','sigfigs=-1')));
end
if ~isempty(h.sform)
  disp(sprintf('sform:'));
  disp(sprintf('%s',mlrnum2str(h.sform,'compact=0','sigfigs=-1')));
end
if ~isempty(h.vol2mag)
  disp(sprintf('vol2mag:'));
  disp(sprintf('%s',mlrnum2str(h.vol2mag,'compact=0','sigfigs=-1')));
end
if ~isempty(h.vol2tal)
  disp(sprintf('vol2tal:'));
  disp(sprintf('%s',mlrnum2str(h.vol2tal,'compact=0','sigfigs=-1')));
end

% get axis information if not passed in
if ~isempty(h.qform)
  axisLabels = mlrImageGetAxisLabels(h.qform);
else
  axisLabels = [];
end

% display axis information
if ~isempty(axisLabels)
  cardinalAxisLabels = {'X','Y','Z'};
  disp(sprintf('Volume orientation is: %s',axisLabels.orient));
  for axisNum = 1:3
    disp(sprintf('Axis %s goes from %s to %s',cardinalAxisLabels{axisNum},axisLabels.dirLabels{axisNum}{1},axisLabels.dirLabels{axisNum}{2}));
  end
end
  
% if there is a talInfo field, display that
if ~isempty(h.talInfo)
  dispHeader(sprintf('%s (Talairach info)',getLastDir(h.filename)));
  disp(sprintf('AC: [%s]',mlrnum2str(h.talInfo.AC,'compact=1','sigfigs=0')));
  disp(sprintf('PC: [%s]',mlrnum2str(h.talInfo.PC,'compact=1','sigfigs=0')));
  disp(sprintf('SAC: [%s]',mlrnum2str(h.talInfo.SAC,'compact=1','sigfigs=0')));
  disp(sprintf('IAC: [%s]',mlrnum2str(h.talInfo.IAC,'compact=1','sigfigs=0')));
  disp(sprintf('PPC: [%s]',mlrnum2str(h.talInfo.PPC,'compact=1','sigfigs=0')));
  disp(sprintf('AAC: [%s]',mlrnum2str(h.talInfo.AAC,'compact=1','sigfigs=0')));
  disp(sprintf('LAC: [%s]',mlrnum2str(h.talInfo.LAC,'compact=1','sigfigs=0')));
  disp(sprintf('RAC: [%s]',mlrnum2str(h.talInfo.RAC,'compact=1','sigfigs=0')));
end


% display detailed header information
if verbose>1
  dispHeader(sprintf('%s (detailed header)',getLastDir(h.filename)));
  hdrFields = fieldnames(h.hdr);
  for iField = 1:length(hdrFields)
    val = h.hdr.(hdrFields{iField});
    if isnumeric(val)
      if (size(val,1) == 1) && (size(val,2) == 1)
	disp(sprintf('%s: %s',hdrFields{iField},mlrnum2str(val)));
      elseif size(val,1) == 1
	disp(sprintf('%s: [%s]',hdrFields{iField},mlrnum2str(val)));
      elseif size(val,2) == 1
	disp(sprintf('%s: [%s]',hdrFields{iField},mlrnum2str(val')));
      else
	disp(sprintf('%s:\n%s',hdrFields{iField},mlrnum2str(val,'compact=0')));
      end
    elseif isstr(val)
      disp(sprintf('%s: %s',hdrFields{iField},val));
    elseif isempty(val)
      disp(sprintf('%s: []',hdrFields{iField}));
    elseif isstruct(val)
      disp(sprintf('%s: struct',hdrFields{iField}));
    else
      disp(sprintf('%s: Unknown type',hdrFields{iField}));
    end      
  end
end

dispHeader;
