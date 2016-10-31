function [tf, view, unknownFields] =  isview(view)
% function [tf view] =  isview(view)
%
% Checks to see if it is a valid view structure. Can be called with
% either one or two output arguments:
%
% tf =  isview(view)
% [tf view] =  isview(view)
%
% tf is logical 1 (true) if view is a valid view structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid view structure by setting optional fields to default
% values.
% 
% djh, 2007

mrGlobals

if (nargout >= 2)
  % Add optional fields and return true if the view with optional fields is
  % valid.
  requiredFields = {'viewNum','viewType','baseVolumes','curBase','curGroup',...
		    'analyses','curAnalysis','ROIs','curROI','prevROIcoords','showROIs','figure','curslice','curScan'};
  optionalFields = {'loadedAnalyses',{};
		    'groupScanNum',[];
		    'labelROIs',0;
		    'roiGroup',{};
		    'sliceOrientation',1;
        'groupSettings',[];
       };
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the view structure it is invalid).
  requiredFields = {'viewNum','viewType','baseVolumes','curBase','curGroup',...
		    'analyses','curAnalysis','ROIs','curROI','prevROIcoords','showROIs','figure','curslice','loadedAnalyses','groupScanNum','labelROIs','curScan','sliceOrientation','roiGroup',...
        'groupSettings'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('view')
  tf = false;
  return
end
if ~isstruct(view)
  tf = false;
  return
end

% Check required fields
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(view,fieldName)
    % mrWarnDlg(['Invalid view, missing field: ',fieldName]);
    tf = false;
  end
end

% Optional fields and defaults
originalView = view;
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(view,fieldName)  
    view.(fieldName) = default;
  end
end
view = orderfields(view);

% Check that viewNum is a number
if ~isfield(view,'viewNum') || ~isnumeric(view.viewNum);
  tf = false;
  return
end

% confirm that there is view in MLR.views with the viewNum
if isempty(view.viewNum) || (view.viewNum < 1) || (view.viewNum > length(MLR.views)) || isempty(MLR.views{view.viewNum})
  tf = false;
  return
end

% Confirm that MLR.views{viewNum} and view have the same fields
names1 = fieldnames(orderfields(MLR.views{view.viewNum}));
names2 = fieldnames(view);
if length(names1) == length(names2)
  tf = all(strcmp(names1,names2));
else
  tf = false;
end

%see if there are any unknown fields
unknownFields = setdiff(fieldnames(originalView),names1);
