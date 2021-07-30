function ea_overlaynii2patch(nii,surfobj,cmap)


% get colors for surface:
bb=nii.mat*[1,size(nii.img,1);1,size(nii.img,2);1,size(nii.img,3);1,1];
[X,Y,Z]=meshgrid(linspace(bb(1,1),bb(1,2),size(nii.img,1)),...
    linspace(bb(2,1),bb(2,2),size(nii.img,2)),...
    linspace(bb(3,1),bb(3,2),size(nii.img,3)));

gradientLevel = length(gray);
defaultColor = [1 1 1]; % Default color for nan values

cmap = [cmap; defaultColor];

ic=isocolors(X,Y,Z,permute(nii.img,[2,1,3]),surfobj.Vertices);
CInd = round(ea_contrast(ic)*gradientLevel+1);
CInd(isnan(CInd)) = gradientLevel + 1; % set to white for now
CData = cmap(CInd,:);

surfobj.FaceVertexCData=CData;
surfobj.FaceColor='interp';
surfobj.EdgeColor='none';
surfobj.SpecularStrength=0.35;
surfobj.SpecularExponent=30;
surfobj.SpecularColorReflectance=0;
surfobj.AmbientStrength=0.07;
surfobj.DiffuseStrength=0.45;


