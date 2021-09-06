function ea_mkdir(pth)
if ~isfolder(pth)
    mkdir(pth);
end
