function vatfv = ea_niiVAT2fvVAT(nii, smooth, smoothingfactor)
if ischar(nii) && isfile(nii)
    nii = ea_load_nii(nii);
end

start = nii.mat * [1;1;1;1];
niisize = size(nii.img);
stop = nii.mat * [niisize,1]';
gv = cell(1,3);
for dim=1:3
    gv{dim}=linspace(start(dim),stop(dim),niisize(dim));
end

[xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});
vatfv=isosurface(xg,yg,zg,permute(nii.img,[2,1,3]),0.5);
caps=isocaps(xg,yg,zg,permute(nii.img,[2,1,3]),0.5);
vatfv.faces=[vatfv.faces;caps.faces+size(vatfv.vertices,1)];
vatfv.vertices=[vatfv.vertices;caps.vertices];

if nargin<2
    smooth = 1;
end

if nargin < 3
    smoothingfactor = 35;
end

if smooth
    try
        vatfv=ea_smoothpatch(vatfv,1,smoothingfactor);
    catch
        try
            cd([ea_getearoot,'ext_libs',filesep,'smoothpatch']);
            mex ea_smoothpatch_curvature_double.c -v
            mex ea_smoothpatch_inversedistance_double.c -v
            mex ea_vertex_neighbours_double.c -v
            vatfv=ea_smoothpatch(vatfv,1,smoothingfactor);
        catch
            warndlg('Patch could not be smoothed. Please supply a compatible Matlab compiler to smooth VTAs.');
        end
    end
end
