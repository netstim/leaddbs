function vatfv = ea_niiVAT2fvVAT(nii)
start = nii.mat * [1;1;1;1];
niisize = size(nii.img);
stop = nii.mat * [niisize,1]';
for dim=1:3
    gv{dim}=linspace(start(dim),stop(dim),niisize(dim));
end
[xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});
vatfv=isosurface(xg,yg,zg,permute(nii.img,[2,1,3]),0.75);
caps=isocaps(xg,yg,zg,permute(nii.img,[2,1,3]),0.5);
vatfv.faces=[vatfv.faces;caps.faces+size(vatfv.vertices,1)];
vatfv.vertices=[vatfv.vertices;caps.vertices];
end