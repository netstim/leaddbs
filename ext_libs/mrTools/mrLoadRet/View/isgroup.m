function [tf group] =  isgroup(group)
% function [tf group] =  isgroup(group)
%
% Checks to see if it is a valid group structure. Can be called with
% either one or two output arguments:
%
% tf =  isgroup(group)
% [tf group] =  isgroup(group)
%
% tf is logical 1 (true) if group is a valid group structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid group structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the group with optional
  % fields is valid.
  requiredFields = {'name','scanParams'};
  optionalFields = {'auxParams',[];};
else
  % Return 0 if the group structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'name','scanParams','auxParams'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('group')
    tf = false;
    return
end
if ~isstruct(group)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(group,fieldName)
		% mrWarnDlg(['Invalid group, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(group,fieldName)  
    group.(fieldName) = default;
  end
end
group = orderfields(group);
