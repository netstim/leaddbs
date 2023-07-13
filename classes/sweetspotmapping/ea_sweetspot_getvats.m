function vatlist = ea_sweetspot_getvats(obj)
% Return list of VATs

numPatient = length(obj.allpatients);
vatlist = cell(numPatient*2,2);

modelLabel = ea_simModel2Label(obj.M.vatmodel);

disp('Construct VAT list...')
for sub=1:numPatient
    % Original VAT E-field
    stimFolder = [obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(0), 'gs_', obj.M.guid];
    try
        vatlist(sub,1) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-R\.nii$'], 0);
    catch
        ea_cprintf('CmdWinWarnings', 'Right side VTA doesn''t exist under stimulation folder:\n%s\n', stimFolder);
        vatlist(sub,1) = {''};
    end
    try
        vatlist(sub,2) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-L\.nii$'], 0);
    catch
        ea_cprintf('CmdWinWarnings', 'Left side VTA doesn''t exist under stimulation folder:\n%s\n', stimFolder);
        vatlist(sub,2) = {''};
    end

    % Mirrored VAT E-field
    ea_genflippedjointnii(vatlist{sub,1}, vatlist{sub,2});
    try
        vatlist(numPatient+sub, 1) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-R_hemidesc-FlippedFromLeft\.nii$'], 0);
    catch
        ea_cprintf('CmdWinWarnings', 'Right side VTA (flipped from left) doesn''t exist under stimulation folder:\n%s\n', stimFolder);
        vatlist(numPatient+sub, 1) = {''};
    end
    try
        vatlist(numPatient+sub, 2) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-L_hemidesc-FlippedFromRight\.nii$'], 0);
    catch
        ea_cprintf('CmdWinWarnings', 'Left side VTA (flipped from right) doesn''t exist under stimulation folder:\n%s\n', stimFolder);
        vatlist(numPatient+sub, 2) = {''};
    end
end
