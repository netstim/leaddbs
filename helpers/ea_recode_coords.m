function ea_recode_coords(options)


[~,~,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
[coords_mm,trajectory,markers]=ea_resolvecoords(markers,options,0);

ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options);

