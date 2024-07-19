function [S] = ea_process_programmer_group(importedPatientS)
    S = importedPatientS.S;
    dt = datetime('now');
    % Convert the datetime object to a string in the format 'yyyymmddHHMMSS'
    formattedDate = datestr(dt, 'yyyymmddHHMMSS');
    S.label = strcat('gs_', formattedDate);
    % S.model = 'OSS-DBS (Butenko 2020)';
    numRows = size(S.activecontacts, 2);
    numCols = size(S.activecontacts, 2);
    
    % Initialize the cell array
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
    % Given data
    S.activecontacts = {secondTerm, firstTerm};  % Assuming these are arrays
    
    % Process each term in the cell array
    for i = 1:length(S.activecontacts)
        term = S.activecontacts{i};  % Extract the term (either secondTerm or firstTerm)
        
        % Check if any element is greater than 1 and change it to 1
        term(term > 1) = 1;
        
        % Save the modified term back to the cell array
        S.activecontacts{i} = term;
    end
    
    % S.label = formattedDate;
    %             S.model = S_old.model;
    S.amplitude = {S.amplitude.rightAmplitude.', S.amplitude.leftAmplitude.'};
    S.monopolarmodel = 0;
    S.sources=[1, 2, 3, 4];
end