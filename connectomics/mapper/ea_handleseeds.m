function [sfile,roilist]=ea_handleseeds(sfile)

% handle seed files
if iscell(sfile) % already supplied in cell format
    if length(sfile)>1
        roilist=1;
    else
        roilist=0;
    end
    
else
    [pth,fn,ext]=fileparts(sfile);
    if strcmp(ext,'.txt')
        roilist=1;
        
        sfile=ea_getrois(sfile);
    else
        roilist=0;
        sfile={sfile};
    end
end
