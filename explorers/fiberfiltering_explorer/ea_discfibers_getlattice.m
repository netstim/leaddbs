function [vatlist,FilesExist] = ea_discfibers_getlattice(obj)
% Return list of VATs

numPatient = length(obj.allpatients);
vatlist = cell(numPatient,2);

modelLabel = ea_simModel2Label(obj.M.vatmodel);

disp('Construct VAT list...')
for sub=1:numPatient
    % Original VAT E-field
    stimFolder = [obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(obj.native), 'gs_', obj.M.guid];
        vatlist{sub,1} = fullfile(stimFolder,'Results_lh','E_field_solution_Lattice.nii');
        FilesExist(sub,1)=exist(vatlist{sub,1},'file');
        vatlist{sub,2} = fullfile(stimFolder,'Results_rh','E_field_solution_Lattice.nii');
        FilesExist(sub,2)=exist(vatlist{sub,2},'file');
end
FilesExist=double(logical(FilesExist)); % convert to zeros and ones.