function ea_dumpnormmethod(options,normmethod,whatmethod)


switch whatmethod
    case 'normmethod'
        varname='norm_method_applied';
    case 'coregctmethod'
        varname='coregct_method_applied';
    case 'coregmrmethod'
        varname='coregmr_method_applied';
end


try load([options.root,options.patientname,filesep,'ea_',whatmethod,'_applied']);
    vn=eval(varname);
catch
    vn={};
end
if ~iscell(vn)
    vc{1}=vn;
    vn=vc;
end

if exist('vn','var')
    try
        vn{end+1}=normmethod;
    catch
        clear(['norm_',whatmethod,'_applied']);
        vn{1}=normmethod;
    end
else
    vn{1}=normmethod;
end
eval([varname,'=vn;']);
save([options.root,options.patientname,filesep,'ea_',whatmethod,'_applied'],varname);
