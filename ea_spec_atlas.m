function ea_spec_atlas(obj, command, cmap, setinterpol)

set(obj, 'LineWidth', 0.1);

if ~setinterpol
    fc=get(obj, 'FaceColor');
    set(obj, 'EdgeColor', fc/1.2);
    colormap(gray);
else
    set(obj, 'FaceColor', 'interp');
    colormap(gray);
end

len = get(obj, 'CData');

switch command
    case 'labeling'
        set(obj, 'FaceVertexAlphaData', repmat(0.5,length(len),1));
        set(obj, 'FaceAlpha', 0.1);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        % set(atls, 'LineStyle', '--');
        % set(atls, 'SpecularColorReflectance', 0);
        % set(atls, 'SpecularExponent', 0);
        % set(atls, 'SpecularStrength', 0)
    case 'mask.nii'
        set(obj, 'FaceVertexAlphaData', repmat(0.5,length(len),1));
        set(obj, 'FaceAlpha', 0.2);
    case 'cortex.nii'
        set(obj, 'FaceVertexAlphaData', repmat(0.5,length(len),1));
        set(obj, 'FaceAlpha', 0.1);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        set(obj, 'LineStyle', '--');
        set(obj, 'SpecularColorReflectance', 0);
        set(obj, 'SpecularExponent', 10);
        set(obj, 'SpecularStrength', 0.5)
    case 'vat'
        try % only works for volumetrics/patches
            set(obj, 'FaceVertexAlphaData', repmat(0.3,length(len),1));
        end
        set(obj, 'FaceAlpha', 0.5);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        set(obj, 'LineStyle', '--');
        set(obj, 'SpecularColorReflectance', 1);
        set(obj, 'SpecularExponent', 3);
        set(obj, 'SpecularStrength', 0.3)
        set(obj, 'DiffuseStrength', 0.4)
        set(obj, 'FaceVertexCData', repmat([0.9,0.1,0.1],size(get(obj,'Vertices'),1),1));
    otherwise
        try % only works for volumetrics/patches
            set(obj, 'FaceVertexAlphaData', repmat(0.3,length(len),1));
        end
        set(obj, 'FaceAlpha', 0.5);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        set(obj, 'LineStyle', '--');
        set(obj, 'SpecularColorReflectance', 1);
        set(obj, 'SpecularExponent', 3);
        set(obj, 'SpecularStrength', 0.3)
        set(obj, 'DiffuseStrength', 0.4)
        % alpha(0.7);
end

%camlight right
