% mlrAnatDBSubjectID.m
%
%        $Id:$ 
%      usage: subjectID = mlrAnatDBSubjectID(v)
%         by: justin gardner
%       date: 06/23/15
%    purpose: Returns a subjectID. This can be used in two ways
%             You can get the subject ID from an MLR session:
%
%             v = newView;
%             subjectID = mlrAnatDBSubjectID(v);
%
%             Or you can get subjectID properly formatted from a number
%             or string
%
%             subjectID = mlrAnatDBSubjectID(25);
%             subjectID = mlrAnatDBSubjectID('s25');
%
function subjectID = mlrAnatDBSubjectID(v)

subjectID = '';
% check if we are passed in a string
if isstr(v)
  subjectID = v;
  % see if it is of the form snnn
  if (length(subjectID) > 1) 
    [idLocStart idLocEnd] = regexp(subjectID,'s\d+');
    id = 0;
    if ~isempty(idLocStart)
      id = str2num(subjectID(idLocStart+1:idLocEnd));
    else
      [idLocStart idLocEnd] = regexp(subjectID,'\d+');
      if ~isempty(idLocStart)
	id = str2num(subjectID(idLocStart:idLocEnd));
      end
    end
    % format subject ID
    subjectID = sprintf('s%04i',id);
  elseif ~isempty(regexp(subjectID,'\d+'))
    subjectID = mlrAnatDBSubjectID(str2num(subjectID));
  end
  return
% check if we are passed in a number
elseif isnumeric(v)
  subjectID = mlrAnatDBSubjectID(sprintf('s%i',v));
  return
% otherwise it should have been a view
elseif isempty(v)
  return
end

if ~isview(v)
  disp(sprintf('(mlrAnatDBSubjectID) Passed in variable is not a view'));
  return
end

% check if we are in the reop
if mlrAnatDBInLocalRepo(v)
  % then just get subject ID from path
  subjectID = getLastDir(fileparts(fileparts(viewGet(v,'homeDir'))));
  subjectID = mlrAnatDBSubjectID(subjectID);
  return
end
    
% get the subject
subjectID = viewGet(v,'subject');
% fromat it
subjectID = mlrAnatDBSubjectID(subjectID);

% ask if this is the correct subjectID
paramsInfo{1} = {'subjectID',subjectID,'The subject ID that this session will be filed under in the anatDB. Usually this is of the form sXXXX. If you do not know it, you may be able to look it up using mglSetSID if you are usuing the mgl ID database'};
params = mrParamsDialog(paramsInfo,'Set subjectID');
if isempty(params)
  subjectID = '';
  return
end
subjectID = params.subjectID;
% check input format
subjectID = mlrAnatDBSubjectID(subjectID);

