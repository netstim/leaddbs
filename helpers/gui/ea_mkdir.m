function ea_mkdir(pth)
if ~isfolder(pth, 'dir')
    mkdir(pth);
end
