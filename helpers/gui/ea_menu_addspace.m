function ea_menu_addspace(handles)

space_menu = uimenu('Label','Space');
change_space_menu = uimenu(space_menu,'Label','Change current space to');

spaces = dir(fullfile(ea_getearoot,'templates','space','*','ea_space_def.mat'));
[~,spaces] = cellfun(@(x) fileparts(strrep(x,'_',' ')), {spaces.folder}', 'uni', false);
sspaces = cellfun(@(x) ea_underscore2space(x), spaces, 'uni', false);

current_space = ea_getspace;
current_space_index = find(cellfun(@(x) strcmp(x, current_space), spaces));

for i = 1:length(spaces)
    c = uimenu(change_space_menu,'Label',sspaces{i},'Callback',{@ea_switchspace,spaces{i}});
    if i == current_space_index
        c.Checked = 'on';
    end
end

uimenu(space_menu,'Label','(Re-)generate aux files for current space, using selected atlas','Callback',{@ea_genauxspace,handles});

import_atlas_menu = uimenu(space_menu,'Label','Import atlases from...');
import_parcellation_menu = uimenu(space_menu,'Label','Import whole-brain parcellations from...');
import_both_menu = uimenu(space_menu,'Label','Import atlases and whole-brain-parcellations from...');
import_custom_menu = uimenu(space_menu,'Label','Import custom .nii file from...');

for i = 1:length(spaces)
    if i == current_space_index
        continue
    end
    uimenu(import_atlas_menu,'Label',sspaces{i},'Callback',{@ea_importspaceassets,spaces{i},'atlases'});
    uimenu(import_parcellation_menu,'Label',sspaces{i},'Callback',{@ea_importspaceassets,spaces{i},'labeling'});
    uimenu(import_both_menu,'Label',sspaces{i},'Callback',{@ea_importspaceassets,spaces{i},'both'});
    uimenu(import_custom_menu,'Label',sspaces{i},'Callback',{@ea_importspaceassets,spaces{i},'custom'});
end

