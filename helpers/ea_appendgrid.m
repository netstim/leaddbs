function  [postops,gfis]=ea_appendgrid(options,postops,gfis,fullpath)
if options.prefs.normalize.createwarpgrids
    if exist([options.root,options.patientname,filesep,'grid.nii'],'file')
        if fullpath
            addstr=[options.root,options.patientname,filesep];
        else
            addstr='';
        end
        postops{end+1}=[addstr,'grid.nii'];
        gfis{end+1}=[addstr,'glgrid.nii'];
%         lfis{end+1}=[addstr,'lgrid.nii'];
    end
end


