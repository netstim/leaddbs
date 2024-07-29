function ea_mirrorfig2blender(fig)

if ~exist('fig','var')
    fig=gcf;
end


% patches
patchObjects = findobj(fig, 'Type', 'patch');
blenderfig = ea_bfigure();

for p=1:length(patchObjects)

% Use the patch method to plot the data in Blender
blenderfig.patch('Vertices', patchObjects(p).Vertices, 'Faces', patchObjects(p).Faces);

end

