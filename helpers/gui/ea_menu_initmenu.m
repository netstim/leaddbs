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
    
    
    
    if ismember('dbs',cmd)
        dbs=uimenu(f,'Label','DBS');
        uimenu(dbs,'Label','Recalculate DBS reconstruction in template space','Callback',{@ea_recalc_reco,handles});
    end

    if ismember('surfice',cmd)
       si=uimenu(f,'Label','Surfice'); 
        uimenu(si,'Label','Visualize DBS-scene in Surfice (template space)','Callback',{@ea_elvis_surfice,handles,0});
        uimenu(si,'Label','Visualize Atlas set in Surfice (template space)','Callback',{@ea_atlvis_surfice,handles,0});
        %uimenu(si,'Label','Visualize DBS-scene in Surfice (native space)','Callback',{@ea_elvis_surfice,handles,1});
        sini=uimenu(si,'Label','Export heatmaps from nifti file(s)');
        siaco=uimenu(sini,'Label','Auto Window');
        uimenu(siaco,'Label','Right hemisphere views','Callback',{@ea_surfice_heatmap_menu,handles,1,0});
        uimenu(siaco,'Label','Left hemisphere views','Callback',{@ea_surfice_heatmap_menu,handles,2,0});
        uimenu(siaco,'Label','Bilateral views','Callback',{@ea_surfice_heatmap_menu,handles,[1,2],0});
        simco=uimenu(sini,'Label','Manual Window');
        uimenu(simco,'Label','Right hemisphere views','Callback',{@ea_surfice_heatmap_menu,handles,1,1});
        uimenu(simco,'Label','Left hemisphere views','Callback',{@ea_surfice_heatmap_menu,handles,2,1});
        uimenu(simco,'Label','Bilateral views','Callback',{@ea_surfice_heatmap_menu,handles,[1,2],1});
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
    
    if ismember('methods',cmd)
       m=uimenu(f,'Label','Methods');
       m_c=uimenu(m,'Label','Show methods popup on completed tasks.','Callback',{@ea_toggle_methods});
       prefs=ea_prefs;
       if prefs.machine.methods.show
           m_c.Checked='on';
       else
           m_c.Checked='off';
       end
       
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
        
       insit(l)=uimenu(g,'Label',[list{l}],'Callback',{@ea_menuinstall,commands{l}});
       if ea_checkinstall(commands{l},1)
           insit(l).Checked='on';
        else
            insit(l).Checked='off';
        end
    end

    % mark that menu has already been installed.
        setappdata(handles.leadfigure,'menuprobe',1);
end


function ea_toggle_methods(hobj,~,~)
switch hobj.Checked
    case 'on'
        ea_setprefs('methods_show',0);
        hobj.Checked='off';
    case 'off'
        ea_setprefs('methods_show',1);
        hobj.Checked='on';
end

