function ea_generate_electrode_specs

[~,matfn]=ea_resolve_elspec;
for elmodel=1:length(matfn)
    try
        feval(['ea_elspec_',matfn{elmodel}],'silent');
    catch
        warning(['Could not generate electrode specification for ',matfn{elmodel},'.']);
    end
end
