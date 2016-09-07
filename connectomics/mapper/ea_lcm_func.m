function ea_lcm_func(options)


if strcmp(options.lcm.func.connectome,'No functional connectome found.')
    return
end
disp('Running functional connectivity...');

if iscell(options.lcm.seeds);
   if length(options.lcm.seeds)==1
      options.lcm.seeds=options.lcm.seeds{1}; % supply as char 
   end
end


cs_fmri_conseed(ea_getconnectomebase,options.lcm.func.connectome,...
    options.lcm.seeds,...
    ea_lcm_resolvecmd(options.lcm.cmd),...
    '0',...
    options.lcm.odir,...
    options.lcm.omask);
disp('Done.');