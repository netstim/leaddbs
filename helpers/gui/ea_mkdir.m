function ea_mkdir(pth)
if ~exist(pth,'dir')
    mkdir(pth);
end