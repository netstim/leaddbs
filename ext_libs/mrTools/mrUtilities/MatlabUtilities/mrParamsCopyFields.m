% mrParamsCopyFields.m
%
%      usage: params2 = mrParamsCopyFields(params1,params2)
%         by: julien besle 
%       date: 04/12/2010
%    purpose: copies (and replaces) fields from structure params1 to structure params2
%              only fields with a corresponding paramInfo cell are copied
%              $Id$

function params2 = mrParamsCopyFields(params1,params2,fieldNames)

if ieNotDefined('params2')
  params2=[];
end

if ~isfield(params1,'paramInfo')
  mrWarnDlg('(mrParamsCopyFields) Structure to copy is not a parameter structure. Use copyFields instead');
  return;
end

if ~ieNotDefined('fieldNames')
  [~,paramsToCopy] = ismember(fieldNames,fieldnames(params1));
  if ~all(paramsToCopy)
    mrWarnDlg(sprintf('(mrParamsCopyFields) Fields %s do not exist',char(fieldNames(~paramsToCopy))'));
    paramsToCopy(~paramsToCopy)=[];
  end
else
  paramsToCopy = 1:length(params1.paramInfo);
end

if ~isfield(params2,'paramInfo')
  params2.paramInfo = [];
  foundParam = 0;
end

for iParam = paramsToCopy
  params2.(params1.paramInfo{iParam}{1})=params1.(params1.paramInfo{iParam}{1});
  %find the paraminfo cell in struct2
  for jParam = 1:length(params2.paramInfo)
    foundParam = 0;
    if strcmp(params2.paramInfo{jParam}{1},params1.paramInfo{iParam}{1})
      params2.paramInfo{jParam}=params1.paramInfo{iParam};
      foundParam=1;
      break;
    end
  end
  if ~foundParam
    params2.paramInfo{end+1} = params1.paramInfo{iParam};
  end
end





