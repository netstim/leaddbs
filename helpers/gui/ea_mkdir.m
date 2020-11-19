function ea_mkdir(pth)
if ~exist(pth,'dir')
    try
        mkdir(pth);
    end
end
