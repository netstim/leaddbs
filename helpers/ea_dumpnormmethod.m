function ea_dumpnormmethod(options,normmethod)


try load([options.root,options.patientname,filesep,'ea_normmethod_applied']); end
if exist('norm_method_applied','var')
    try
        norm_method_applied{end+1}=normmethod;
    catch
        clear norm_method_applied
        norm_method_applied{1}=normmethod;
    end
else
    norm_method_applied{1}=normmethod;
end
save([options.root,options.patientname,filesep,'ea_normmethod_applied'],'norm_method_applied');
