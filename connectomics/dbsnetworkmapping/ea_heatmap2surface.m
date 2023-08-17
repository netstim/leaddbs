function h=ea_heatmap2surface(nii,surface,sides,colormap,opts)

% surface can be a char (smoothed/full) or a cell pointing to rh and lh
% .mz3 surface files.
if ~exist('opts','var')
    opts.posvisible=1;
    opts.negvisible=1;
end

if ~exist('surface','var')
    surface='smoothed';
end
if ~exist('sides','var')
    sides=1:2;
end
if ischar(surface)
    % first draw correct surface
    switch lower(surface)
        case 'smoothed'
            if ismember(1,sides); [rh.faces, rh.vertices] = ea_readMz3([ea_space,'surf_smoothed.rh.mz3']); end
            if ismember(2,sides); [lh.faces, lh.vertices] = ea_readMz3([ea_space,'surf_smoothed.lh.mz3']); end
        case 'full'
            if ismember(1,sides); [rh.faces, rh.vertices] = ea_readMz3([ea_space,'surf.rh.mz3']); end
            if ismember(2,sides); [lh.faces, lh.vertices] = ea_readMz3([ea_space,'surf.lh.mz3']); end
    end
elseif    iscell(surface) % assuming rh and lh
    if ismember(1,sides); [rh.faces, rh.vertices] = ea_readMz3(surface{1}); end
    if ismember(2,sides); [lh.faces, lh.vertices] = ea_readMz3(surface{2}); end
    surface='custom';
end

if ~exist('colormap','var')
    colormap=ea_redblue;
end
gradientLevel = length(colormap);

% Check cmap
defaultColor = [1 1 1]; % Default color for nan values
cmap = [colormap; defaultColor];

if ischar(nii)
    res=ea_load_nii(nii);
else
    res=nii;
end

% get colors for surface:
bb=res.mat*[1,size(res.img,1);1,size(res.img,2);1,size(res.img,3);1,1];
[X,Y,Z]=meshgrid(linspace(bb(1,1),bb(1,2),size(res.img,1)),...
    linspace(bb(2,1),bb(2,2),size(res.img,2)),...
    linspace(bb(3,1),bb(3,2),size(res.img,3)));

if ~opts.posvisible
    res.img(res.img>0)=0;
end

if ~opts.negvisible
    res.img(res.img<0)=0;
end

[pth,niiname]=fileparts(res.fname);
h = {};
if ismember(1,sides)
    ic=isocolors(X,Y,Z,permute(res.img,[2,1,3]),rh.vertices);
    if any(~isnan(ic))
        CInd = round(ea_contrast(ic)*gradientLevel+1);
        CInd(isnan(CInd)) = gradientLevel + 1; % set to white for now
        rhCData = cmap(CInd,:);
    else
        rhCData = repmat(defaultColor, size(rh.vertices,1), 1);
    end
    h{1}=patch('Faces',rh.faces,'Vertices',rh.vertices,'FaceColor','interp','EdgeColor','none','FaceVertexCData',rhCData,...
        'SpecularStrength',0.35,'SpecularExponent',30,'SpecularColorReflectance',0,'AmbientStrength',0.07,'DiffuseStrength',0.45,'FaceLighting','gouraud');
    h{1}.Tag=[niiname,'_right'];
end

if ismember(2,sides)
    ic=isocolors(X,Y,Z,permute(res.img,[2,1,3]),lh.vertices);
    if any(~isnan(ic))
        CInd = round(ea_contrast(ic)*gradientLevel+1);
        CInd(isnan(CInd)) = gradientLevel + 1; % set to white for now
        lhCData = cmap(CInd,:);
    else
        lhCData = repmat(defaultColor, size(lh.vertices,1), 1);
    end
    h{2}=patch('Faces',lh.faces,'Vertices',lh.vertices,'FaceColor','interp','EdgeColor','none','FaceVertexCData',lhCData,...
        'SpecularStrength',0.35,'SpecularExponent',30,'SpecularColorReflectance',0,'AmbientStrength',0.07,'DiffuseStrength',0.45,'FaceLighting','gouraud');
    h{2}.Tag=[niiname,'_left'];
end

