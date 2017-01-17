function ea_menu_addspace(handles)

g = uimenu('Label','Space');
s = uimenu(g,'Label','Change current space to');
spaces=dir([ea_getearoot,'templates',filesep,'space',filesep]);
for space=1:length(spaces)
    if ~strcmp(spaces(space).name(1),'.')
        spacename=spaces(space).name;
        if strcmp(spacename,ea_getspace)
            spacename=['-> ',spacename];
        end
        sspacename=ea_sub2space(spacename);
        c=uimenu(s,'Label',sspacename,'Callback',{@ea_switchspace,spacename});
    end
end
s = uimenu(g,'Label','(Re-)generate aux files for current space, using selected atlas','Callback',{@ea_genauxspace,spacename});