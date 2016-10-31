function [tf roi] =  isroi(roi)
% function [tf roi] =  isroi(roi)
%
% Checks to see if it is a valid roi structure. Can be called with
% either one or two output arguments:
%
% tf =  isroi(roi)
% [tf roi] =  isroi(roi)
%
% tf is logical 1 (true) if roi is a valid roi structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid roi structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the roi with optional
  % fields is valid.
  requiredFields = {'name','voxelSize','xform'};
  optionalFields = {'viewType',[];
		    'coords',[];
		    'date',datestr(now);
		    'color','blue';
		    'notes','';
		    'sformCode',1;
		    'vol2mag',[];
		    'vol2tal',[];
		    'createdBy','';
		    'createdOnBase','';
		    'displayOnBase','';
		    'createdFromSession','';
		    'branchNum',[];
		    'subjectID',[];
		   };
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'color','coords','date','name','viewType','voxelSize','xform','sformCode','vol2mag','vol2tal','createdBy','createdOnBase','displayOnBase','createdFromSession','branchNum','subjectID'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('roi')
    tf = false;
    return
end
if ~isstruct(roi)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(roi,fieldName)
		% mrWarnDlg(['Invalid roi, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(roi,fieldName)  
    roi.(fieldName) = default;
  end
end

% remove any fields that are not required or optional
if nargout == 2
  roiFieldNames = fieldnames(roi);
  for f = 1:length(roiFieldNames)
    % check required fields
    if ~any(strcmp(roiFieldNames{f},requiredFields))
      % check optional fields, (only check first column of
      % field names, not the default values...
      match = strcmp(roiFieldNames{f},optionalFields);
      if ~any(match(:,1))
	roi = rmfield(roi,roiFieldNames{f});
      end
    end
  end
end


% if roi coords doesn't have the fourth dimension then
% set it to 1
if tf && ~isempty(roi.coords) && (size(roi.coords,1) == 3)
  roi.coords(4,:) = 1;
end

% make sure voxelSize is
if isfield(roi,'voxelSize')
  roi.voxelSize = roi.voxelSize(:)';
end
% order the fields
roi = orderfields(roi);
