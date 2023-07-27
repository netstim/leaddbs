function filepath = ea_path_helper(filepath)
% handle special characters in the path for cli compatibility

if ischar(filepath) % support char or cell input
    filepath = {filepath};
    waschar = 1;
else
    waschar = 0;
end

for i=1:length(filepath)
    filepath{i} = GetFullPath(filepath{i});

    if ispc
        filepath{i} = ['"',filepath{i},'"'];
    else
        filepath{i} = regexprep(filepath{i},'[[](){}''"; &]', '\\$0');
    end
end

if waschar
    filepath = filepath{1};
end
