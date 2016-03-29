function ea_spec_atlas(atls,command,cmap,setinterpol)

% switch command
%     
%     case 'gp.nii'
%         if ~setinterpol
%             set(atls, 'FaceColor', [0.6,0.6,0])
%         end
%         
%     case 'stn.nii'
%         if ~setinterpol
%             set(atls, 'FaceColor', colorc);
%         end
%     case 'thalamus.nii'
%         if ~setinterpol
%             set(atls, 'FaceColor', colorc);
%         end
%         
%     case 'putamen.nii'
%         if ~setinterpol
%             set(atls, 'FaceColor', colorc);
%         end
%         
%     case 'labeling'
%         if ~setinterpol
%             set(atls, 'FaceColor', colorc);
%         end
%         
% end
set(atls, 'LineWidth',0.1);

if ~setinterpol
    fc=get(atls,'FaceColor');
    set(atls, 'EdgeColor',fc/1.2);
    colormap(cmap)
else
    set(atls,'FaceColor','interp');
    colormap(cmap)
end
len=get(atls,'CData');


switch command
    case 'labeling'
set(atls,'FaceVertexAlphaData',repmat(0.5,length(len),1));
set(atls,'FaceAlpha',0.1);


set(atls, 'EdgeColor','none');

set(atls, 'FaceLighting', 'gouraud');
%set(atls, 'LineStyle', '--');
%set(atls, 'SpecularColorReflectance', 0);
%set(atls, 'SpecularExponent', 0);
%set(atls, 'SpecularStrength', 0)
case 'mask.nii'
set(atls,'FaceVertexAlphaData',repmat(0.5,length(len),1));
set(atls,'FaceAlpha',0.2);
    case 'cortex.nii'
set(atls,'FaceVertexAlphaData',repmat(0.5,length(len),1));
set(atls,'FaceAlpha',0.1);

set(atls, 'EdgeColor','none');

set(atls, 'FaceLighting', 'gouraud');
set(atls, 'LineStyle', '--');
set(atls, 'SpecularColorReflectance', 0);
set(atls, 'SpecularExponent', 10);
set(atls, 'SpecularStrength', 0.5)
    otherwise   
set(atls,'FaceVertexAlphaData',repmat(0.3,length(len),1));
set(atls,'FaceAlpha',0.5);



set(atls, 'EdgeColor','none');

set(atls, 'FaceLighting', 'gouraud');
set(atls, 'LineStyle', '--');
set(atls, 'SpecularColorReflectance', 0);
set(atls, 'SpecularExponent', 10);
set(atls, 'SpecularStrength', 0.5)
%alpha(0.7);
end




%camlight right
