function h=ea_parcellation2surface(nii,surface,sides,colormap)

% surface can be a char (smoothed/full) or a cell pointing to rh and lh
% .mz3 surface files.


if ~exist('surface','var')
    surface='auto';
end
if ~exist('sides','var')
    sides=1:2;
end

if ~exist('colormap','var')
    colormap=ea_color_wes('all');
end
%gradientLevel = length(colormap);

expansionfactor=1; % extrapolate to reduce artifacts
expcolormap=zeros(size(colormap,1)*expansionfactor,3);
count=1;
for c=1:expansionfactor:length(expcolormap)
expcolormap(c:c+expansionfactor-1,:)=repmat(colormap(count,:),expansionfactor,1);
count=count+1;
end





if ischar(nii)
    if exist(nii,'file')
    res=ea_load_nii(nii);
    else
        try
            res=ea_load_nii(ea_niigz(fullfile(ea_space([],'labeling'),[nii,'.nii'])));
        catch
            ea_error('Input File does not exist.');
        end
    end
else
    res=nii;
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
        case 'auto'

            resp=permute(res.img,[2,1,3]);

            resp=resp>0;
            resp=smooth3(double(resp),'gaussian',[5,5,5]);
            resp=resp>0.95;
            [x,y,z,D] = subvolume(resp,[1,size(resp,2),1,round(size(resp,1)/2),1,size(resp,3)]);
            lh=isosurface(x,y,z,D, 0.5);
            [x,y,z,D] = subvolume(resp,[1,size(resp,2),round(size(resp,1)/2),size(resp,1),1,size(resp,3)]);
            rh=isosurface(x,y,z,D, 0.5);
          

            lh.vertices=res.mat*[lh.vertices';ones(1,size(lh.vertices,1))];
            lh.vertices=lh.vertices(1:3,:)';
            lh=ea_smoothpatch(lh,1,3);
            rh.vertices=res.mat*[rh.vertices';ones(1,size(rh.vertices,1))];
            rh.vertices=rh.vertices(1:3,:)';
            rh=ea_smoothpatch(rh,1,3);
    end
elseif    iscell(surface) % assuming rh and lh
    if ismember(1,sides); [rh.faces, rh.vertices] = ea_readMz3(surface{1}); end
    if ismember(2,sides); [lh.faces, lh.vertices] = ea_readMz3(surface{2}); end
    surface='custom';
end







res.img=round(res.img).*expansionfactor;
while max(res.img(:))>size(expcolormap,1)
    expcolormap=[expcolormap;expcolormap];
end

% Check cmap
defaultColor = [1 1 1]; % Default color for nan values
cmap = [expcolormap; defaultColor];

% get colors for surface:
bb=res.mat*[1,size(res.img,1);1,size(res.img,2);1,size(res.img,3);1,1];
[X,Y,Z]=meshgrid(linspace(bb(1,1),bb(1,2),size(res.img,1)),...
    linspace(bb(2,1),bb(2,2),size(res.img,2)),...
    linspace(bb(3,1),bb(3,2),size(res.img,3)));


[pth,niiname]=fileparts(res.fname);
h = {};
if ismember(1,sides)
    ic=isocolors(X,Y,Z,permute(res.img,[2,1,3]),rh.vertices);
    if any(~isnan(ic))
        CInd = round((ic)+1);
        CInd(isnan(CInd)) = size(cmap,1); % set to white for now
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
        CInd = round((ic)+1);
        CInd(isnan(CInd)) = size(cmap,1); % set to white for now
        lhCData = cmap(CInd,:);
    else
        lhCData = repmat(defaultColor, size(lh.vertices,1), 1);
    end
    h{2}=patch('Faces',lh.faces,'Vertices',lh.vertices,'FaceColor','interp','EdgeColor','none','FaceVertexCData',lhCData,...
        'SpecularStrength',0.35,'SpecularExponent',30,'SpecularColorReflectance',0,'AmbientStrength',0.07,'DiffuseStrength',0.45,'FaceLighting','gouraud');
    h{2}.Tag=[niiname,'_left'];
end

