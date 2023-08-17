function ea_mkdir(pth)

if ~iscell(pth)
    pth = {pth};
end

cellfun(@(x) ~isfolder(x) && mkdir(x), pth);
