function outputfolder = ea_getoutputfolder(sfile,con)
if strcmp(con, 'wFTR.mat')
    con = 'Patient-specific fiber tracts';
else
    con = strrep(con, '>', '_');
end

outputfolder = [fileparts(sfile{1}), filesep, con, filesep];
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder);
end
