function pth=ea_path_formatlab(pth)
%this makes the path safe to handle for matlab tools and derivatives (which
%would fail otherwise (e.g. mkdir, ea_mkdir)
%this method is needed on windows platform (ispc) if ea_path_helper() is used

%stripping string wrapped with quotes ""
pth = ea_stripquotes(pth);
