function ea_menu_addvats(toolsmenu, handles)

v = uimenu(toolsmenu,'Label','VAT');
uimenu(v,'Label','Aggregate all files of selected stimulation (Left VAT)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_hemi-L')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right VAT)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_hemi-R')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left E-Field)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-efield_hemi-L')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right E-Field)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-efield_hemi-R')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left normalized E-Field)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-efieldgauss_hemi-R')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right normalized E-Field)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-efieldgauss_hemi-R')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left fMRI seed)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_seed-fMRI_hemi-L')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right fMRI seed)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_seed-fMRI_hemi-R')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Combined fMRI seed)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_seed-fMRI')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left dMRI seed)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_seed-dMRI_hemi-L')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right dMRI seed)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_seed-dMRI_hemi-R')});
uimenu(v,'Label','Aggregate all files of selected stimulation (Combined dMRI seed)','Callback',{@(src, evt) ea_aggregateVTA(handles,'sim-binary_seed-dMRI')});
