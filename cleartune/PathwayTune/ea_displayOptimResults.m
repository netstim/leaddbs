function ea_displayOptimResults(activ_profile_fig, symptom_profile_fig)
    %info1 = [ea_getearoot, 'ext_libs/OSS-DBS/OSSSettings_pt1.png'];
    %info2 = [ea_getearoot, 'ext_libs/OSS-DBS/OSSSettings_pt2.png'];
    if ispc
        winopen(activ_profile_fig);
        winopen(symptom_profile_fig);
    elseif ismac
        system(['open ', activ_profile_fig]);
        system(['open ', symptom_profile_fig]);
    else
        if ~system('which xdg-open >/dev/null 2>&1')
            system(['xdg-open ', symptom_profile_fig]);
            system(['xdg-open ', activ_profile_fig]);
        elseif ~system('which display >/dev/null 2>&1')
            system(['display -resize 35% ', symptom_profile_fig, ' &']);
            system(['display -resize 35% ', activ_profile_fig, ' &']);
        else
            figure('NumberTitle', 'off', 'MenuBar', 'none');
            imshow(info2);
            figure('NumberTitle', 'off', 'MenuBar', 'none');
            imshow(info1);
        end
    end
    
end