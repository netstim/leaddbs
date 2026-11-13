function ea_spec_atlas(obj, command, cmap, setinterpol, uniformColor)
% The optional uniformColor skips per-vertex color/alpha data for a patch
% whose displayed color is constant. The caller sets the scene colormap.
uniform = nargin >= 5 && ~isempty(uniformColor);

set(obj, 'LineWidth', 0.1);
if uniform
    set(obj, 'FaceColor', uniformColor);
end

if ~setinterpol
    fc=get(obj, 'FaceColor');
    set(obj, 'EdgeColor', fc/1.2);
    if ~uniform
        colormap(gray);
    end
else
    if ~uniform
        set(obj, 'FaceColor', 'interp');
        colormap(gray);
    end
end

if ~uniform
    len = get(obj, 'CData');
end

switch command
    case 'labeling'
        if ~uniform
            set(obj, 'FaceVertexAlphaData', repmat(0.5,length(len),1));
        end
        set(obj, 'FaceAlpha', 0.1);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        % set(atls, 'LineStyle', '--');
        % set(atls, 'SpecularColorReflectance', 0);
        % set(atls, 'SpecularExponent', 0);
        % set(atls, 'SpecularStrength', 0)
    case 'mask.nii'
        if ~uniform
            set(obj, 'FaceVertexAlphaData', repmat(0.5,length(len),1));
        end
        set(obj, 'FaceAlpha', 0.2);
    case 'cortex.nii'
        if ~uniform
            set(obj, 'FaceVertexAlphaData', repmat(0.5,length(len),1));
        end
        set(obj, 'FaceAlpha', 0.1);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        set(obj, 'LineStyle', '--');
        set(obj, 'SpecularColorReflectance', 0);
        set(obj, 'SpecularExponent', 10);
        set(obj, 'SpecularStrength', 0.5)
    case 'vat'
        try % only works for volumetrics/patches
            if ~uniform
                set(obj, 'FaceVertexAlphaData', repmat(0.3,length(len),1));
            end
        end
        set(obj, 'FaceAlpha', 0.5);
        set(obj, 'EdgeColor', 'none');
        set(obj, 'FaceLighting', 'gouraud');
        set(obj, 'LineStyle', '--');
        set(obj, 'SpecularColorReflectance', 1);
        set(obj, 'SpecularExponent', 3);
        set(obj, 'SpecularStrength', 0.3)
        set(obj, 'DiffuseStrength', 0.4)
        if uniform
            set(obj, 'FaceColor', [0.9,0.1,0.1]);
        else
            set(obj, 'FaceVertexCData', repmat([0.9,0.1,0.1],size(get(obj,'Vertices'),1),1));
        end
    otherwise
        try % only works for volumetrics/patches
            if ~uniform
                set(obj, 'FaceVertexAlphaData', repmat(0.3,length(len),1));
            end
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
