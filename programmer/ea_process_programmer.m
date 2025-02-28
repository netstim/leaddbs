function [S] = ea_process_programmer(options)
    %% Initialize Paths, Read in data from programmer through data.json, Delete helper files (inputData, data)

    stimDir = fullfile(options.subj.stimDir, ea_getspace);

    importedS = loadjson(fullfile(stimDir, 'data.json'));
    
    ea_delete(fullfile(stimDir, 'data.json'));
    ea_delete(fullfile(stimDir, 'inputData.json'));

    if isfield(importedS, 'message')
        S = importedS;
        return;
    end

    %% Reorganize input for programmer to support lead-dbs convention

    S = importedS.S;
    numRows = size(S.activecontacts, 2);
    numCols = size(S.activecontacts, 2);
    
    newVariable = cell(2, 4);
    
    % Fill the cell array
    for i = 1:2
        for j = 1:4
            newVariable{i, j} = S.activecontacts((i-1)*4+j, :);
        end
    end

    S.activecontacts = newVariable;

    firstTerm=S.activecontacts{1,1} + S.activecontacts{1,2} + S.activecontacts{1,3} + S.activecontacts{1,4};
    secondTerm=S.activecontacts{2,1} + S.activecontacts{2,2} + S.activecontacts{2,3} + S.activecontacts{2,4};

    S.activecontacts={secondTerm, firstTerm};
    S.activecontacts = {secondTerm, firstTerm};
    
    for i = 1:length(S.activecontacts)
        term = S.activecontacts{i};
        
        term(term > 1) = 1;
        
        S.activecontacts{i} = term;
    end

    S.amplitude = {S.amplitude.rightAmplitude.', S.amplitude.leftAmplitude.'};
    S.monopolarmodel = 0;
    S.sources=[1, 2, 3, 4];
end