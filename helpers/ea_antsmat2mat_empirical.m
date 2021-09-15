function emat=ea_antsmat2mat_empirical(options)
% Extra empirical check of the scrf matrix

moving = options.subj.brainshift.anat.moving;
anchor = options.subj.brainshift.anat.anchor;

mov = ea_load_nii(moving);
[vxx, vyy, vzz] = ind2sub(size(mov.img),1:100:numel(mov.img));
XYZ_mov_vx = [vxx; vyy; vzz; ones(1,length(vxx))];
XYZ_fix_mm = ea_map_coords(XYZ_mov_vx(1:3,:), moving, options.subj.brainshift.transform.instore, anchor, 'ANTS');

voxmov2mmfix = (XYZ_mov_vx(1:4,:)'\[XYZ_fix_mm;ones(1,size(XYZ_fix_mm,2))]')';

emat = voxmov2mmfix/mov.mat; % combine mats
