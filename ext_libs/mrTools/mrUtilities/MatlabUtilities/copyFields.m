% copyFields.m
%
%      usage: struct2 = copyFields(struct1,struct2,indices)
%         by: julien besle 
%       date: 03/12/2010
%    purpose: copies (and replaces) all fields from struct1 to struct2, or struct2(indices) if indices is defined
%              $Id$

function struct2 = copyFields(struct1,struct2, indices)


if ieNotDefined('struct2')
  struct2=[];
end
if ieNotDefined('indices')
  indices=1;
end

fieldNames=fieldnames(struct1);

%if struct2 is an array of structures, first check that all fields exist in struct2
if length(struct2)>1
  struct2FieldNames=fieldnames(struct1);
  fieldsToCreate = find(~ismember(fieldNames,struct2FieldNames));
  if ~isempty(fieldsToCreate)
    for iField = 1:length(fieldNames)
      struct2 = set(struct2,fieldsToCreate,fieldNames{fieldsToCreate},{1},[]);
    end
  end
end
    
 
for iStruct = indices
  for iField = 1:length(fieldNames)
    struct2(iStruct).(fieldNames{iField})=struct1.(fieldNames{iField});
  end
end
