function ea_menu_initmenu(handles,cmd)

callingfunction=dbstack;
callingfunction=callingfunction(4).name;
menuprobe=getappdata(handles.leadfigure,'menuprobe');
if isempty(menuprobe)
    % tools menu  & edit prefs present in all apps.
    f = uimenu('Label','Tools');
    uimenu(f,'Label','Edit Preferences File...','Callback',{@ea_editprefs});
    uimenu(f,'Label','Reset Preferences to Default...','Callback',{@ea_restoreprefs});
    if ismember('acpc',cmd)
        uimenu(f,'Label','Convert ACPC/MNI coordinates (Horn 2017)','Callback',{@ea_acpcquery,handles.leadfigure});
    end

    if ismember('export',cmd)
        e = uimenu(f,'Label','Export');
        uimenu(e,'Label','Export .PDF files for selected patient(s)','Callback',{@ea_exportpat,'PDF',handles});
        uimenu(e,'Label','Export .STL files for selected patient(s)','Callback',{@ea_exportpat,'STL',handles});
        uimenu(e,'Label','Export .PLY files for selected patient(s)','Callback',{@ea_exportpat,'PLY',handles});

        d = uimenu(f,'Label','Convert');
        uimenu(d,'Label','Convert selected atlas to .STL','Callback',{@ea_exportatlas,'STL',handles});
        uimenu(d,'Label','Convert selected atlas to .PLY','Callback',{@ea_exportatlas,'PLY',handles});
    end

    if ismember('applynorm',cmd)
        uimenu(f,'Label','Apply Patient Normalization to file...','Callback',{@ea_applynormtofile_menu,handles,0});
        uimenu(f,'Label','Apply Patient Inverse Normalization to file...','Callback',{@ea_applynormtofile_menu,handles,1});

    end
    if ismember('cluster',cmd)
        ea_menu_addsubmit(handles);
    end

    if ismember('space',cmd)
        ea_menu_addspace(handles);
    end

    if ismember('checkregfigs',cmd)
       uimenu(f,'Label','Generate Checkreg figures','Callback',{@ea_gencheckreg,handles});
    end

    if ismember('transfer',cmd)
       ea_menu_addtransfer(handles,callingfunction);
    end

    if ismember('vats',cmd)
       ea_menu_addvats(f,handles);
    end

    % always add install addons
    g = uimenu('Label','Install');
    [list,commands]=ea_checkinstall('list');
    for l=1:length(list)
        if ea_checkinstall(commands{l},1)
           addstr='v';
        else
            addstr='ø';
        end
       uimenu(g,'Label',[list{l},' (',addstr,')'],'Callback',{@ea_menuinstall,commands{l}});
    end

    % mark that menu has already been installed.
        setappdata(handles.leadfigure,'menuprobe',1);
end
