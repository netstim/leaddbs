function emat=ea_antsmat2mat_empirical(directory)

if directory(end) ~= filesep
    directory = [directory, filesep];
end
movim=[directory,'scrf',filesep,'movim.nii'];
[options.root,options.patientname]=fileparts(directory);
options.root = [options.root, filesep];
options.prefs=ea_prefs;
[~,anatpresent]=ea_assignpretra(options);

fixim=[directory,'scrf',filesep,anatpresent{1}];

mov=ea_load_nii(movim);
[vxx,vyy,vzz]=ind2sub(size(mov.img),1:100:numel(mov.img));
XYZ_mov_vx=[vxx;vyy;vzz;ones(1,length(vxx))];
XYZ_fix_mm = ea_map_coords(XYZ_mov_vx(1:3,:), movim, [directory,'scrf',filesep,'scrf_instore.mat'], fixim, 'ANTS');

voxmov2mmfix=(XYZ_mov_vx(1:4,:)'\[XYZ_fix_mm;ones(1,size(XYZ_fix_mm,2))]')';

emat=voxmov2mmfix/mov.mat; % combine mats


