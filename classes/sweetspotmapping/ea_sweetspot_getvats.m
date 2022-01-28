function vatlist = ea_sweetspot_getvats(obj)
% Return list of VATs

if obj.M.ui.detached
    pthprefix = [fileparts(obj.leadgroup),filesep];
else
    pthprefix = '';
end

numPatient = length(obj.allpatients);
vatlist = cell(numPatient*2,2);

modelLabel = ea_simModel2Label(obj.M.vatmodel);

disp('Construct VAT list...')
for sub=1:numPatient
    % Original VAT E-field
    stimFolder = [pthprefix, obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(0), 'gs_', obj.M.guid];
    vatlist(sub,1) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-R\.nii$'], 0);
    vatlist(sub,2) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-L\.nii$'], 0);

    % Mirrored VAT E-field
    ea_genflippedjointnii(vatlist{sub,1}, vatlist{sub,2});
    vatlist(numPatient+sub, 1) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-R_hemidesc-FlippedFromLeft\.nii$'], 0);
    vatlist(numPatient+sub, 2) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-L_hemidesc-FlippedFromRight\.nii$'], 0);
end
