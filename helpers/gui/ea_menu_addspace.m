function ea_menu_addspace(handles)

g = uimenu('Label','Space');
s = uimenu(g,'Label','Change current space to');
spaces=dir([ea_getearoot,'templates',filesep,'space',filesep]);
for space=1:length(spaces)
    if ~strcmp(spaces(space).name(1),'.')
        spacename=spaces(space).name;
        sspacename=ea_underscore2space(spacename);
        c=uimenu(s,'Label',sspacename,'Callback',{@ea_switchspace,spacename});
        if strcmp(spacename,ea_getspace)
            c.Checked='on';
        end
    end
end
j = uimenu(g,'Label','(Re-)generate aux files for current space, using selected atlas','Callback',{@ea_genauxspace,handles});

k1 = uimenu(g,'Label','Import atlases from...');
for space=1:length(spaces)
    if ~strcmp(spaces(space).name(1),'.')
        spacename=spaces(space).name;
        if strcmp(spacename,ea_getspace)
            continue
        end
        sspacename=ea_underscore2space(spacename);
        c=uimenu(k1,'Label',sspacename,'Callback',{@ea_importspaceassets,spacename,'atlases'});
    end
end

k2 = uimenu(g,'Label','Import whole-brain parcellations from...');
for space=1:length(spaces)
    if ~strcmp(spaces(space).name(1),'.')
        spacename=spaces(space).name;
        if strcmp(spacename,ea_getspace)
            continue
        end
        sspacename=ea_underscore2space(spacename);
        c=uimenu(k2,'Label',sspacename,'Callback',{@ea_importspaceassets,spacename,'labeling'});
    end
end

k3 = uimenu(g,'Label','Import atlases and whole-brain-parcellations from...');
for space=1:length(spaces)
    if ~strcmp(spaces(space).name(1),'.')
        spacename=spaces(space).name;
        if strcmp(spacename,ea_getspace)
            continue
        end
        sspacename=ea_underscore2space(spacename);
        c=uimenu(k3,'Label',sspacename,'Callback',{@ea_importspaceassets,spacename,'both'});
    end
end

k3 = uimenu(g,'Label','Import custom .nii file from...');
for space=1:length(spaces)
    if ~strcmp(spaces(space).name(1),'.')
        spacename=spaces(space).name;
        if strcmp(spacename,ea_getspace)
            continue
        end
        sspacename=ea_underscore2space(spacename);
        c=uimenu(k3,'Label',sspacename,'Callback',{@ea_importspaceassets,spacename,'custom'});
    end
end
