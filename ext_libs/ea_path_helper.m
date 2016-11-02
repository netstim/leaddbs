function filepath = ea_path_helper(filepath)
% handle special characters in the path for cli compatibility

if isempty(fileparts(filepath)) && ~strcmp(filepath,'.')
    filepath = ['.', filesep, filepath];
end

if ispc
    filepath = ['"',filepath,'"'];
else
    filepath = regexprep(filepath,'[() &]', '\\$0');
end
