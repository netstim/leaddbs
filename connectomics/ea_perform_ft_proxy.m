function ea_perform_ft_proxy(options)
% function ea_ft_mesotracking_reisert
% function ea_ft_deterministictracking_kroon
% function ea_ft_globaltracking_reisert
% function ea_ft_gqi_yeh

eval([options.lc.struc.ft.method,'(options)']); % triggers the fibertracking function and passes the options struct to it.

% BIDS: Save method log in connectomics/log/
directory=[options.root,options.patientname,filesep];
logDir = fullfile(directory, 'connectomics', 'log');
if ~exist(logDir, 'dir')
    mkdir(logDir);
end

try load(fullfile(logDir, 'ea_ftmethod_applied.mat')); end
if exist('ft_method_applied','var')
    try
        ft_method_applied{end+1}=options.lc.struc.ft.method;
    catch
        clear ft_method_applied
        ft_method_applied{1}=options.lc.struc.ft.method;
    end
else
    ft_method_applied{1}=options.lc.struc.ft.method;
end
ft_method_applied=options.lc.struc.ft.method;

save(fullfile(logDir, 'ea_ftmethod_applied.mat'),'ft_method_applied');
