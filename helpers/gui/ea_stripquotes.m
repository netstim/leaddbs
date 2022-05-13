function pth=ea_stripquotes(pth)
%stripping string wrapped with quotes ""

pth = strip(pth,'both','"');
