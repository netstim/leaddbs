function spec_atlas(atls,colorc,command,setinterpol)

switch command
    
    case 'gp.nii'
        if ~setinterpol
            set(atls, 'FaceColor', [0.6,0.6,0])
        end
        
    case 'stn.nii'
        if ~setinterpol
            set(atls, 'FaceColor', colorc);
        end
    case 'thalamus.nii'
        if ~setinterpol
            set(atls, 'FaceColor', colorc);
        end
        
    case 'putamen.nii'
        if ~setinterpol
            set(atls, 'FaceColor', colorc);
        end
        
    case 'aal'
        if ~setinterpol
            set(atls, 'FaceColor', colorc);
        end
        
end
set(atls, 'LineWidth',0.1);

if ~setinterpol
    fc=get(atls,'FaceColor');
    set(atls, 'EdgeColor',fc/1.2);
    colormap jet
else
    set(atls,'FaceColor','interp');
    colormap jet
end
len=get(atls,'CData');


switch command
    case 'aal'
set(atls,'FaceVertexAlphaData',repmat(0.5,length(len),1));
set(atls,'FaceAlpha',0.05);
    otherwise   
set(atls,'FaceVertexAlphaData',repmat(0.3,length(len),1));
set(atls,'FaceAlpha',0.7);
%alpha(0.7);
end





set(atls, 'FaceLighting', 'gouraud');
set(atls, 'LineStyle', '--');
set(atls, 'SpecularColorReflectance', 0);
set(atls, 'SpecularExponent', 10);
set(atls, 'SpecularStrength', 1)

%camlight right
shading interp