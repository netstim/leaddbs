
function out = dti_realigndef_ui(P)

mat = load(P.matname{1},'M');

fi = realign_def(P.yname{1},P.iyname{1},inv(mat.M));

out.finame_y = fi.y;
out.finame_iy = fi.iy;


function out = realign_def(yfi,iyfi,mat)


[yDef,ymat] = get_def(yfi);
[iyDef,iymat] = get_def(iyfi);

% realign y field
ymat = mat*ymat;

% realign iy field
niyDef{1} = single(iyDef{1}*mat(1,1) + iyDef{2}*mat(1,2) + iyDef{3}*mat(1,3) + mat(1,4));
niyDef{2} = single(iyDef{1}*mat(2,1) + iyDef{2}*mat(2,2) + iyDef{3}*mat(2,3) + mat(2,4));
niyDef{3} = single(iyDef{1}*mat(3,1) + iyDef{2}*mat(3,2) + iyDef{3}*mat(3,3) + mat(3,4));
iyDef = niyDef;


% save y field
[p n e] = fileparts(yfi);
out.y = save_def(yDef,ymat,[ n],p);

% save iy field
[p n e] = fileparts(iyfi);
out.iy = save_def(iyDef,iymat,[ n],p);








function [Def,mat] = get_def(file)
% Load a deformation field saved as an image

P      = [repmat(file,3,1), [',1,1';',1,2';',1,3']];
V      = spm_vol(P);
Def    = cell(3,1);
Def{1} = spm_load_float(V(1));
Def{2} = spm_load_float(V(2));
Def{3} = spm_load_float(V(3));
mat    = V(1).mat;






function fname = save_def(Def,mat,ofname,odir)
% Save a deformation field as an image

if isempty(ofname), fname = {}; return; end;

fname = {fullfile(odir,[ofname '.nii'])};
dim   = [size(Def{1},1) size(Def{1},2) size(Def{1},3) 1 3];
dtype = 'FLOAT32';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(fname{1},dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'VECTOR';
N.intent.name = 'Mapping';
N.descrip = 'Deformation field';
create(N);
N.dat(:,:,:,1,1) = Def{1};
N.dat(:,:,:,1,2) = Def{2};
N.dat(:,:,:,1,3) = Def{3};
return;