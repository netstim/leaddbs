function thresh=ea_detthresh(atlases,atlas,img)
% function that will return the threshold according to user settings for a
% specific atlas.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if isfield(atlases,'threshold')
    switch atlases.threshold.type
        case 'percentage'
            try
                sso=sort(img(img>0));
                thresh=sso(round(length(sso)*(1-atlases.threshold.value))); % preserve % of voxels.
            catch
                thresh=sso(1);
            end
        case 'percentage_vector'
            sso=sort(img(img>0));
            try
                thresh=sso(round(length(sso)*(1-atlases.threshold.value(atlas)))); % preserve % of voxels.
            catch
                thresh=sso(1);
            end
        case 'relative_intensity'
            
            thresh=max(img(:))*(1-atlases.threshold.value);
        case 'relative_intensity_vector'
            thresh=max(img(:))*(1-atlases.threshold.value(atlas));
        case 'absolute_intensity'
            thresh=atlases.threshold.value;
        case 'absolute_intensity_vector'
            thresh=atlases.threshold.value(atlas);
        otherwise
            warning(['Threshold type not recognized: ',atlases.threshold.type,'. Overwriting with default.']);
            thresh=max(img(:))*0.5;
            
    end
    
else
    
    thresh=max(img(:))*0.5;
end