function ea_visprogrammer(resultfig, options, S, elstruct)

    funcs = ea_regexpdir(ea_getearoot, 'ea_genvat_.*\.m$', 0);
    funcs = regexp(funcs, '(ea_genvat_.*)(?=\.m)', 'match', 'once');
    names = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);
    
    setappdata(resultfig,'genvatfunctions',funcs);
    setappdata(resultfig,'vatfunctionnames',names);
    vfnames=getappdata(resultfig,'vatfunctionnames');
    
    [~,ix]=ismember(S.model,vfnames);
    vfs=getappdata(resultfig,'genvatfunctions');
    try
        ea_genvat=eval(['@',vfs{ix}]);
    catch
        keyboard
    end
    if isequal(S.model, 'OSS-DBS (Butenko 2020)') 
            [~, stimparams] = ea_genvat_butenko(S, options, resultfig);
            flix = 1; 
    else
        for side=1:2
            try 
    %             [vtafv, vtavolume] = ea_genvat_horn(elstruct.coords_mm, S, side, options, S.label, resultfig);
    %             [vtafv,vtavolume] = feval(ea_genvat,coords,M.S(pt),side,options,['gs_',M.guid],handles.leadfigure);
                [vtafv,vtavolume] = feval(ea_genvat,elstruct.coords_mm,S,side,options,S.label,resultfig);
                vtaCalcPassed(side) = 1;
            catch 
                vtafv=[];
                vtavolume=0;
                vatCalcPassed(side) = 0;
            end
            stimparams(1,side).VAT(1).VAT = vtafv;
            stimparams(1,side).volume = vtavolume;
        end
    end           
    
    setappdata(resultfig,'stimparams',stimparams);
    setappdata(resultfig,'curS',S);
    hmchanged = 1;
    ea_calc_vatstats(resultfig,options,hmchanged);
    input_file_path = strcat(options.earoot, 'programmer/inputData.json');
    fid = fopen(input_file_path, 'w');
    fclose(fid);

end