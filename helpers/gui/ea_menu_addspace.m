function ea_menu_addspace(handles)

s = uimenu('Label','Space');
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