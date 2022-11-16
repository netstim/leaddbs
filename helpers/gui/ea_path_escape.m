function [pth]=ea_path_escape(pth)
%escaping string path, by making special characters printable by printf format
% warning #1: pay attention to not pass a string that does NOT need to escape special characters.
% warning #2: make sure to NOT double escape a string

    if iscell(pth)
        %assume cell array of string paths
        pth = cellfun(@escape_path_str, pth, 'UniformOutput', false);
    else
        %assume single path
        pth = escape_path_str(pth);
    end
end

function strpth=escape_path_str(strpth)
    strpth = strrep(strpth, '\', '\\'); % backslash
    % strpth = strrep(strpth, '%', '%%'); % percent (percent is not used in a path, DO NOTHING for this)
end