function [groupdir, analysisFile] = ea_genDatasetFromGroupAnalysis(analysisFile)
% Function to create dataset folder from group analysis file (single file
% outside of a dataset)

analysisFile = GetFullPath(analysisFile);

dataset = regexp(analysisFile, '(?<=dataset-)(.+)(?=_analysis-.+\.mat$)', 'match', 'once');
analysis = regexp(analysisFile, '(?<=dataset-.+_analysis-)(.+)(?=\.mat$)', 'match', 'once');
if isfolder(fullfile(fileparts(analysisFile), dataset))
    dataset = inputdlg(sprintf('Folder ''%s'' already exists.\nPlease input a new dataset name:', dataset), 'New Dataset Name', [1 35], {[dataset, '1']});
    if isempty(dataset)
        error('Please input a new dataset name!');
    else
        dataset = dataset{1};
        movefile(analysisFile, regexprep(analysisFile, '(?<=dataset-)(.+)(?=_analysis-.+\.mat$)', dataset));
        analysisFile = regexprep(analysisFile, '(?<=dataset-)(.+)(?=_analysis-.+\.mat$)', dataset);
    end
end

ea_cprintf('CmdWinWarnings', 'Creating new dataset folder: %s\n', fullfile(fileparts(analysisFile), dataset));
groupdir = fullfile(fileparts(analysisFile), dataset, 'derivatives', 'leadgroup', analysis, filesep);

leaddbsFolder = fullfile(fileparts(analysisFile), dataset, 'derivatives', 'leaddbs');
ea_mkdir(groupdir);
ea_mkdir(leaddbsFolder);
fclose(fopen(fullfile(leaddbsFolder, 'Miniset_flag.json'), 'w'));

load(analysisFile, 'M');
M.ui.groupdir = groupdir;
M.root = M.ui.groupdir;
for p=1:length(M.patient.list)
    [oldPatientFolder, patientTag] = fileparts(M.patient.list{p});
    M.patient.list{p} = strrep(M.patient.list{p}, oldPatientFolder, leaddbsFolder);
    if isfield(M, 'stats')
        ea_mkdir(fullfile(leaddbsFolder, patientTag));
        ea_stats = M.stats(p).ea_stats;
        save(fullfile(leaddbsFolder, patientTag, [patientTag, '_desc-stats.mat']), 'ea_stats');
    end

    if isfield(M, 'elstruct')
        ea_mkdir(fullfile(leaddbsFolder, patientTag, 'reconstruction'));
        for e=1:length(M.elstruct(p).coords_mm)
            reco.props(e).elmodel = M.elstruct(p).elmodel;
            reco.props(e).manually_corrected = 1;
            reco.mni.coords_mm(e) = M.elstruct(p).coords_mm(e);
            reco.mni.markers(e) = M.elstruct(p).markers(e);
            reco.mni.trajectory(e) = M.elstruct(p).trajectory(e);
            reco.electrode(e).dbs.elmodel = M.elstruct(p).elmodel;
        end
        save(fullfile(leaddbsFolder, patientTag, 'reconstruction', [patientTag, '_desc-reconstruction.mat']), 'reco');
    end
end
movefile(analysisFile, groupdir);
analysisFile = ea_getGroupAnalysisFile(groupdir);
save(analysisFile, 'M');
