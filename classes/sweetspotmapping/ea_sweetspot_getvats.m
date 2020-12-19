function vatlist = ea_sweetspot_getvats(obj)
% Return list of VATs

if obj.M.ui.detached
    pthprefix = [fileparts(obj.leadgroup),filesep];
else
    pthprefix = '';
end

numPatient = length(obj.allpatients);
vatlist = cell(numPatient*2,2);

disp('Construct VAT list...')
for sub=1:numPatient % Original VAT E-field
    vatlist{sub,1} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'vat_efield_right.nii'];
    vatlist{sub,2} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'vat_efield_left.nii'];
end

for sub=1:numPatient % Mirrored VAT E-field
    ea_genflippedjointnii([pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0) ,'gs_',obj.M.guid,filesep, 'vat_efield_right.nii'],...
        [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'vat_efield_left.nii']);
    vatlist{numPatient+sub,1} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'fl_vat_efield_left.nii'];
    vatlist{numPatient+sub,2} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'fl_vat_efield_right.nii'];
end
