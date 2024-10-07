function axonStateFolder = ea_sourceIndex4AxonStates(outputPaths, side, source_index)

% Add source index to each Axon_state file
% IMPORTANT:this is a special case for multisource prob. PAM
% By Butenko, konstantinmgtu@gmail.com

arguments
    outputPaths         % various paths to conform with lead-dbs BIDS structure 
    side                {mustBeNumeric} % hemisphere index (0 - rh, 1 - lh)
    source_index        % index of the stim. source (not contact!)
end

switch side
    case 0
        sideCode = 'rh';
    case 1
        sideCode = 'lh';
end

% run one source and check AxonStates
axon_state_files = dir([outputPaths.HemiSimFolder,filesep,'Results',filesep,'Axon_state_*']); 
axonStateFolder = [outputPaths.outputDir,filesep,'AxonStates'];
% this might not work for unilateral stim!
if (source_index == 1 || source_index == 5) && side == 0 && isfolder(axonStateFolder)
    ea_delete(axonStateFolder)
    mkdir(axonStateFolder)
elseif ~isfolder(axonStateFolder)
    mkdir(axonStateFolder)
end
for file_i = 1:length(axon_state_files)
    oldfile = [axon_state_files(file_i).folder,filesep,axon_state_files(file_i).name];
    newfile = [axonStateFolder,filesep,axon_state_files(file_i).name(1:end-4),'_',sideCode,num2str(source_index),'.mat'];
    movefile(oldfile,newfile)
end
