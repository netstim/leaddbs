function outputPaths = ea_get_oss_outputPaths(options,S)
% Prepare various paths to conform with Lead-DBS BIDS structure. 
% by Butenko and Li, konstantinmgtu@gmail.com

arguments
    options % lead-dbs options for electrode reconstruction and stimulation
    S       % lead-dbs stimulation settings
end

outputPaths.subDescPrefix = ['sub-', options.subj.subjId, '_desc-'];
subSimPrefix = ['sub-', options.subj.subjId, '_sim-'];
outputPaths.outputDir = [options.subj.stimDir, filesep, ea_nt(options.native), S.label];
outputPaths.outputBasePath = [outputPaths.outputDir, filesep, subSimPrefix];
ea_mkdir(outputPaths.outputDir);
if options.native
    outputPaths.templateOutputDir = [options.subj.stimDir, filesep, ea_nt(0), S.label];
    ea_mkdir(outputPaths.templateOutputDir);
    outputPaths.templateOutputBasePath = [outputPaths.templateOutputDir, filesep, subSimPrefix];
end