function ea_aggregateVTA(handles, exportType)

% Override export type using prefs when exporting dMRI/fMRI seed
if contains(exportType, 'sim-binary_seed-')
    prefs = ea_prefs;
    exportType = strrep(exportType, 'binary', prefs.lcm.vatseed);
end

if ~startsWith(handles.seeddefpopup.String{handles.seeddefpopup.Value}, 'Use VAT: ')
    ea_error('Please select a stimulation first.');
else
    stimLabel = erase(handles.seeddefpopup.String{handles.seeddefpopup.Value}, 'Use VAT: ');
end

uipatdir = getappdata(handles.leadfigure,'uipatdir');
exportDir = uigetdir('','Select where to save files...');

txt = fopen(fullfile(exportDir,'export.txt'), 'w');
for pt=1:length(uipatdir)
    [~, subPrefix] = fileparts([uipatdir{pt}, '_']);

    stimParams = ea_regexpdir(fullfile(uipatdir{pt}, 'stimulations', ea_nt(0), stimLabel), 'stimparameters\.mat$', 0);
    load(stimParams{1}, 'S');
    modelLabel = ea_simModel2Label(S.model);
    exportType = regexprep(exportType, 'sim-([a-z]+)', ['sim-$1_model-', modelLabel]);

    copyfile(fullfile(uipatdir{pt}, 'stimulations', ea_nt(0), stimLabel, [subPrefix, exportType, '.nii']), exportDir);
    fprintf(txt, '%s\n', fullfile(exportDir, [subPrefix, exportType, '.nii']));
end
fclose(txt);

msgbox('Files successfully exported.');
