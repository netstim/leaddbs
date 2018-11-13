function ea_repair_lg_recon(directory)
[coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(directory);
ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,ea_getptopts(directory));