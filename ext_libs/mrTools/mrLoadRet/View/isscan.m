function [tf scanParams] =  isscan(scanParams)
% function [tf scanParams] =  isscan(scanParams)
%
% Checks to see if it is a valid scanParams structure. Can be called with
% either one or two output arguments:
%
% tf =  isscan(scanParams)
% [tf scanParams] =  isscan(scanParams)
%
% tf is logical 1 (true) if scanParams is a valid scanParams structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid scanParams structure by setting optional fields to default
% values.
% 
% djh, 2007
% jlg, 4/2007 check for unknown fields
% djh, 7/2007 allow optional fields with defaults

if (nargout == 2)
  % Add optional fields and return true if the scanParams with optional fields is
  % valid.
  requiredFields = {'description','fileName','fileType','niftiHdr',...
    'voxelSize','totalFrames','junkFrames','nFrames',...
    'dataSize','framePeriod'};
  optionalFields = {'originalGroupName',[];
		    'originalFileName',scanParams.fileName;
		    'totalJunkedFrames',[];
		    'vol2tal',[];
		    'vol2mag',[]};

else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'description','fileName','fileType','niftiHdr',...
    'voxelSize','totalFrames','junkFrames','nFrames',...
    'dataSize','framePeriod','originalFileName','originalGroupName','totalJunkedFrames','vol2tal','vol2mag'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('scanParams')
    tf = false;
    return
end
if ~isstruct(scanParams)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(scanParams,fieldName)
		 mrWarnDlg(['Invalid scanParams, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(scanParams,fieldName)  
    % set argument
    scanParams.(fieldName) = default;
%    mrWarnDlg(sprintf('(isscan) Setting default field %s',fieldName));
  end
end
scanParams = orderfields(scanParams);

% check for unknown fields
allFields = {requiredFields{:}, optionalFields{:}};
scanFields = fieldnames(scanParams);
for f = 1:length(scanFields)
  fieldName = scanFields{f};
  if ~any(strcmp(fieldName,allFields))
%    mrWarnDlg(sprintf('(isscan) Unknown field %s removed',fieldName));
    tf = false;
  end
end
