function fv=ea_fem_getmask(options,nosmooth)

if ~exist('nosmooth','var')
   nosmooth=0;
end
% load nifti
try
    switch options.native
        case 0 % template space
            nii=ea_load_nii(ea_niigz([ea_space(options,'atlases'),options.prefs.machine.vatsettings.horn_atlasset,filesep,'gm_mask.nii.gz']));
        case 1 % native space
            toptions=options;
            toptions.atlasset=options.prefs.machine.vatsettings.horn_atlasset;
            ea_ptspecific_atl(toptions); % make sure atlas has been warped.
            nii=ea_load_nii(ea_niigz([options.root,options.patientname,filesep,'atlases',filesep,options.prefs.machine.vatsettings.horn_atlasset,filesep,'gm_mask.nii.gz']));
    end
catch ME
    ea_cprintf('CmdWinErrors', '\n%s\n\n', ME.message);
    ea_error('Failed to load gray matter mask.', simpleStack = 1);
end

[xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>0)); %(mean(nii.img(nii.img~=0))/3))); % find 3D-points that have correct value.
if isempty(xx)
    fv.vertices=[]; % no gray matter.
    fv.faces=[];
    return
end
XYZ = [xx,yy,zz]; % concatenate points to one matrix.
XYZ = ea_vox2mm(XYZ,nii.mat);

bb=[1,1,1;size(nii.img)];
bb=ea_vox2mm(bb,nii.mat);
gv=cell(3,1);
for dim=1:3
    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(nii.img,dim));
end
[X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
if options.prefs.hullsmooth
    nii.img = smooth3(nii.img,'gaussian',options.prefs.hullsmooth);
end
fv=isosurface(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fvc=isocaps(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
fv.vertices=[fv.vertices;fvc.vertices];

if ~nosmooth
    fv=ea_smoothpatch(fv,[],ceil(options.prefs.hullsmooth/2));
end

% figure
% patch('vertices',fv.vertices,'faces',fv.faces,'facecolor','r');
