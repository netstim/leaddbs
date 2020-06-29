function outputfolder = ea_getoutputfolder(sfile,con)
prefs = ea_prefs;
if strcmp(con, prefs.FTR_normalized) || strcmp(con, strrep(prefs.FTR_unnormalized, '.mat', '_anat.mat'))
    con = 'Patient''s fiber tracts';
else
    con = strrep(con, '>', '_');
end

outputfolder = [fileparts(sfile{1}), filesep, con, filesep];
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder);
end
