function [sfile,roilist]=ea_handleseeds(sfile)

% handle seed files
if iscell(sfile) % already supplied in cell format
    if length(sfile)>1
        roilist=1;
    else
        [~,~,ext]=fileparts(sfile{1});
        if strcmp(ext,'.txt')
            roilist=1;
            sfile=ea_getrois(sfile{1});
        else
            roilist=0;
            sfile=sfile';
        end
    end
else
    [~,~,ext]=fileparts(sfile);
    if strcmp(ext,'.txt')
        roilist=1;
        sfile=ea_getrois(sfile);
    else
        roilist=0;
        sfile={sfile};
    end
end

if size(sfile,1)<size(sfile,2)
    sfile=sfile';
end

for s=1:length(sfile)
   if strcmp(sfile{s,1}(1),'[')
       comma=strfind(sfile{s},',');
       lh=sfile{s}(2:comma-1);
       rh=sfile{s}(comma+1:end-1);
       sfile{s,1}=lh; sfile{s,2}=rh;
   end
end
