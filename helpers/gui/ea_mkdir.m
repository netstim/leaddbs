function ea_mkdir(pth)

%this makes the path safe to handle for matlab (just in case a quoted path
%is passed as input)
pth=ea_path_formatlab(pth);

if ~isfolder(pth)
    mkdir(pth);
end
