function outputfolder=ea_getoutputfolder(sfile,con)
file=sfile{1};
[pth,fn,ext]=fileparts(file); % exit to same folder as seed.
outputfolder=[pth,filesep,con,filesep];
if ~exist(outputfolder,'dir')
    mkdir(outputfolder);
end
