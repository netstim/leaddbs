function filepath = ea_path_helper(filepath)
% handle special characters in the path for cli compatibility

if ischar(filepath) % support char or cell input
    filepath = {filepath};
    waschar = 1;
else
    waschar = 0;
end

for i=1:length(filepath)
    if isempty(fileparts(filepath{i})) && ~strcmp(filepath{i},'.')
        filepath{i} = ['.', filesep, filepath{i}];
    end

    if ispc
        %ensure quotes are not added two times
        if ea_string_startendwithquotes(filepath{i})
            %strip quotes (in the case there was a badly formatted path,
            %with only the quote at the beginning of the string-path
            %(Example: pth='"D:/path/to/dir')
            filepath{i}=ea_stripquotes(filepath{i});
        end
        %add quotes (to take care of spaces in the path)
        filepath{i} = ['"',filepath{i},'"'];
    else
        filepath{i} = regexprep(filepath{i},'[[''() &]]', '\\$0');
    end
end

if waschar
    filepath = filepath{1};
end
