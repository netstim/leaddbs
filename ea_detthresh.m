function thresh=ea_detthresh(atlases,atlas,img)
% Function that will return the threshold according to user settings for a
% specific region in the atlas.
%
% USAGE:
%
%    thresh = ea_detthresh(atlases,atlas,img)
%
% INPUTS:
%    atlases:       structure loaded from 'atlas_index.mat' under atlas folder
%    atlas:         index of query region (needed when the threshold preset is a vector)
%    img:           3D image data of the query region, used to determine the real 
%                   threshold in case of relative threshold used
%
% OUTPUT:
%    thresh:        threshold to be applied to the region defined in the atlas
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Daniel Duarte, Documentation

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
            if atlases.threshold.value(atlas)==1 % include everything
                thresh=sso(end)-eps;
            else
                try
                    thresh=sso(round(length(sso)*(1-atlases.threshold.value(atlas)))); % preserve % of voxels.
                catch
                    thresh=sso(1);
                end
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
