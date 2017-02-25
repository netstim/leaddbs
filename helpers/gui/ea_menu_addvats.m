function ea_menu_addvats(toolsmenu,handles)


v = uimenu(toolsmenu,'Label','VAT');
uimenu(v,'Label','Aggregate all files of selected stimulation (Left VAT)','Callback',{@ea_aggregateVTA,handles,'vat_left'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right VAT)','Callback',{@ea_aggregateVTA,handles,'vat_right'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left E-Field)','Callback',{@ea_aggregateVTA,handles,'vat_efield_left'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right E-Field)','Callback',{@ea_aggregateVTA,handles,'vat_efield_right'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left normalized E-Field)','Callback',{@ea_aggregateVTA,handles,'vat_efield_gauss_left'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right normalized E-Field)','Callback',{@ea_aggregateVTA,handles,'vat_efield_gauss_right'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left fMRI seed)','Callback',{@ea_aggregateVTA,handles,'vat_seed_compound_fMRI_l'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right fMRI seed)','Callback',{@ea_aggregateVTA,handles,'vat_seed_compound_fMRI_r'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Combined fMRI seed)','Callback',{@ea_aggregateVTA,handles,'vat_seed_compound_fMRI'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Left dMRI seed)','Callback',{@ea_aggregateVTA,handles,'vat_seed_compound_dMRI_l'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Right dMRI seed)','Callback',{@ea_aggregateVTA,handles,'vat_seed_compound_dMRI_r'});
uimenu(v,'Label','Aggregate all files of selected stimulation (Combined dMRI seed)','Callback',{@ea_aggregateVTA,handles,'vat_seed_compound_dMRI'});


        
     

