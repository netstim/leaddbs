function filepath = ea_path_helper(filepath)
% handle special characters in the path for cli compatibility

if ispc
    filepath = ['"',filepath,'"'];
else
    filepath = regexprep(filepath,'[() &]', '\\$0');
end

if isempty(fileparts(filepath))
    filepath = ['.', filesep, filepath];
end
