function pth=ea_path_formatlab(pth)
%this makes the path safe to handle for matlab tools and derivatives (which
%would fail otherwise (e.g. mkdir, ea_mkdir)
%this method is needed on windows platform (ispc) if ea_path_helper() is used

    if iscell(pth)
        %assume cell array of string paths
        for str_i=1:length(pth)
            pth{str_i}=fix_path_str(pth{str_i});
        end
    else
        %assume single path
        pth=fix_path_str(pth);
    end
end

function strpth=fix_path_str(strpth)
    %stripping string wrapped with quotes ""
    strpth = ea_stripquotes(strpth);
end

