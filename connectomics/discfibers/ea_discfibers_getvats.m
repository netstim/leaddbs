function [vats, vatlist] = ea_discfibers_getvats(obj)
if exist(fullfile(fileparts(obj.leadgroup), 'disctracts', 'VATs.mat'),'file')
    data = load(fullfile(fileparts(obj.leadgroup), 'disctracts', 'VATs.mat'));
    vats = data.vats;
    vatlist = data.vatlist;
    return;
end

if obj.M.ui.detached
    pthprefix = [fileparts(obj.leadgroup),filesep];
else
    pthprefix = '';
end

numPatient = length(obj.allpatients);
vatlist = cell(numPatient*2,2);

for sub=1:numPatient % all patients - for connected fibers selection ? and always flip
    vatlist{sub,1} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'vat_efield_right.nii'];
    vatlist{sub,2} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'vat_efield_left.nii'];
end

for sub=1:numPatient % all patients - for connected fibers selection ? and always flip
    ea_genflippedjointnii([pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0) ,'gs_',obj.M.guid,filesep, 'vat_efield_right.nii'],...
        [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'vat_efield_left.nii']);
    vatlist{numPatient+sub,1} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'fl_vat_efield_left.nii'];
    vatlist{numPatient+sub,2} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep, 'fl_vat_efield_right.nii'];
end

vats = cellfun(@load_vat, vatlist, 'Uni', 0);

ea_mkdir(fullfile(fileparts(obj.leadgroup),'disctracts'));
save(fullfile(fileparts(obj.leadgroup),'disctracts','VATs.mat'),'vats', 'vatlist', '-v7.3');


function vat = load_vat(vatfile)
% Return a vat struct with file name, affine matrix, dimension and image

vat.fname = vatfile;
vatnii = ea_load_nii(vatfile);
vat.mat = vatnii.mat;
vat.dim = vatnii.dim;
vat.img = vatnii.img;
