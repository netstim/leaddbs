function ea_generate_electrode_specs

[~,matfn]=ea_resolve_elspec;
for elmodel=1:length(matfn)
   feval(['ea_elspec_',matfn{elmodel}],'silent');
end
