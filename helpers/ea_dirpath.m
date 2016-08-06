function allfnames=ea_dirpath(pattern)

pth=path;

[~,seps]=ismember(pth,pathsep);
seps=find(seps);
seps=[0,seps];
for di=1:length(seps)
    try
   fnames=dir([pth(seps(di)+1:seps(di+1)-1),filesep,pattern]);
    catch
           fnames=dir([pth(seps(di)+1:end-1),filesep,pattern]);
    end
    if ~isempty(fnames)
        if ~exist('allfnames','var')
            allfnames=fnames;
        else
            allfnames=[allfnames;fnames];
        end
    end
end
if ~exist('allfnames','var') % return empty dir output.
    allfnames=fnames;
end


