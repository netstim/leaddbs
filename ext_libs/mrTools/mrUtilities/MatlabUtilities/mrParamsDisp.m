% mrParamsDisp.m
%
%      usage: mrParamsDisp(params)
%        $Id$
%         by: julien besle 
%       date: 17/12/2010
%    purpose: displays a params structure using mrParamsDialog. Deals with params substructures
%

function mrParamsDisp(params,title)

if ieNotDefined('title')
  title='Parameters';
end

fieldNames = fieldnames(params);
processed = zeros(length(fieldNames),1);
cParam = 0;

if isfield(params,'paramInfo')
  processed(ismember(fieldNames,'paramInfo'))=1;
  paramInfo = params.paramInfo;
  %find fields that are in paraminfo and that are not empty
  for iParam = 1:length(paramInfo)
    if isfield(params,paramInfo{iParam}{1})
      processed(ismember(fieldNames,paramInfo{iParam}{1}))=1;
      if ~isempty(params.(paramInfo{iParam}{1}))
        cParam = cParam+1;
        newParamInfo{cParam} = cell(1,4);
        %find the help info
        helpInfo = '';
        for j=3:length(paramInfo{iParam})
          if isempty(paramInfo{iParam}{j}) || ((length(strfind(paramInfo{iParam}{j},'=')) ~= 1) && (~isempty(strfind(paramInfo{iParam}{j},' ')))) || ...
                ~isempty(strfind(paramInfo{iParam}{j}(1:strfind(paramInfo{iParam}{j},'=')),' '))
            helpInfo = paramInfo{iParam}{j};
          end
        end
        %just keep the first field (name of parameter)
        newParamInfo{cParam} = paramInfo{iParam}(1);
        %copy the value
        newParamInfo{cParam}{2} = params.(paramInfo{iParam}{1});
        %the third field is the help
        newParamInfo{cParam}{3} = helpInfo;
        %make it non-editable
        newParamInfo{cParam}{4} = 'editable=0';
      end
    end
  end
end

%find additional substructures with paraminfo fields
subParams = [];
for iField = 1:length(fieldNames)
  if ~processed(iField)
    %transform field into cell to simplify handling cellarray of structure
    if ~iscell(params.(fieldNames{iField}))
      thisField={params.(fieldNames{iField})};
    else
      thisField=params.(fieldNames{iField});
    end
    for iCell = 1:length(thisField)
      if isstruct(thisField{iCell}) 
        %if at least one cell is a structure, then we consider it a parameter substructure
        processed(iField)=1;
        if length(thisField)==1
          thisTitle = fieldNames{iField};
        else
          thisTitle = sprintf('%s{%i}',fieldNames{iField},iCell);
        end
        subParams{end+1} = {thisTitle,0,'type=pushbutton',sprintf('buttonString=Show %s',thisTitle),'callback',{@mrParamsDisp,thisField{iCell},thisTitle}};
      end
    end
  end
end

%find additional non-empty parameters
for iField = 1:length(fieldNames)
  if ~processed(iField) && ~isempty(params.(fieldNames{iField}))
    cParam = cParam+1;
    newParamInfo{cParam} = cell(1,4);
    %just keep the first field (name of parameter)
    newParamInfo{cParam}{1} = fieldNames{iField};
    %copy the value
    newParamInfo{cParam}{2} = params.(fieldNames{iField});
    %the third field is the help
    newParamInfo{cParam}{3} = 'No help available';
    %make it non-editable
    newParamInfo{cParam}{4} = 'editable=0';
  end
end

for iParam = 1:length(newParamInfo)
    %and the second field is the actual value
  if islogical(newParamInfo{iParam}{2})
    if newParamInfo{iParam}{2}
      newParamInfo{iParam}{2} = 'True';
    else
      newParamInfo{iParam}{2} = 'False';
    end
  elseif iscell(newParamInfo{iParam}{2})
    newParamInfo{iParam}{5} = 'type=stringArray';
  end
end

newParamInfo = [newParamInfo subParams];

mrParamsDialog(newParamInfo,title,'modal=0');

%this is a version where parameters are only processed if they also are is the paramInfo field
% 
% if isfield(params,'paramInfo')
%   paramInfo = params.paramInfo;
%   cParam = 0;
%   %find fields that are in paraminfo and that are not empty
%   for iParam = 1:length(paramInfo)
%     if ~isfield(params,paramInfo{iParam}{1})
%       mrWarnDlg(sprintf('(mrParamsDisp) field %s is undefined',paramInfo{iParam}{1}));
%     else
%       if ~isempty(params.(paramInfo{iParam}{1}))
%         cParam = cParam+1;
%         newParamInfo{cParam} = cell(1,4);
%         %find the help info
%         helpInfo = '';
%         for j=3:length(paramInfo{iParam})
%           if isempty(paramInfo{iParam}{j}) || ((length(strfind(paramInfo{iParam}{j},'=')) ~= 1) && (~isempty(strfind(paramInfo{iParam}{j},' ')))) || ...
%                 ~isempty(strfind(paramInfo{iParam}{j}(1:strfind(paramInfo{iParam}{j},'=')),' '))
%             helpInfo = paramInfo{iParam}{j};
%           end
%         end
%         %just keep the first field
%         newParamInfo{cParam} = paramInfo{iParam}(1);
%         %the third field is the help
%         newParamInfo{cParam}{3} = helpInfo;
%         %make it non-editable
%         newParamInfo{cParam}{4} = 'editable=0';
%         %and the seond field is the actual value
%         if islogical(params.(paramInfo{iParam}{1}))
%           if params.(paramInfo{iParam}{1})
%             newParamInfo{cParam}{2} = 'True';
%           else
%             newParamInfo{cParam}{2} = 'False';
%           end
%         elseif iscell(params.(paramInfo{iParam}{1}))
%           newParamInfo{cParam}{2} = params.(paramInfo{iParam}{1});
%           newParamInfo{cParam}{5} = 'type=stringArray';
%         else
%           newParamInfo{cParam}{2} = params.(paramInfo{iParam}{1});
%         end
%       end
%     end
%   end
%   %find additional substructures with paraminfo fields
%   fieldNames = fieldnames(params);
%   for iField = 1:length(fieldNames)
%     %transform field into cell to simplify handling cellarray of structu
%     if ~iscell(params.(fieldNames{iField}))
%       params.(fieldNames{iField})={params.(fieldNames{iField})};
%     end
%     for iCell = 1:length(params.(fieldNames{iField}))
%       if isstruct(params.(fieldNames{iField}){iCell}) 
%         if length(params.(fieldNames{iField}))==1
%           thisTitle = fieldNames{iField};
%         else
%           thisTitle = sprintf('%s{%i}',fieldNames{iField},iCell);
%         end
%         if isfield(params.(fieldNames{iField}){iCell},'paramInfo')
%           newParamInfo{end+1} = {thisTitle,0,'type=pushbutton',sprintf('buttonString=Show %s',thisTitle),'callback',{@mrParamsDisp,params.(fieldNames{iField}){iCell},thisTitle}};
%         else
%           newParamInfo{end+1} = {thisTitle,'This is not a parameter structure','type=statictext'};
%         end
%       end
%     end
%   end
% 
%   mrParamsDialog(newParamInfo,title,'modal=0');
% else
%   mrWarnDlg('(mrParamsDisp) This is not a parameter structure');
% end
% 
