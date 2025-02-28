function pos_and_neg = ea_check_for_SSE_tracts(connectomeFolder,slider_names)

% check for which symptoms negative tracts were also exported

pos_and_neg = zeros(size(slider_names,2),1);
pathways_files = dir([connectomeFolder,filesep,'*.mat']);
for slider_i = 1:size(slider_names,2)
    for pathway_i = 1:size(pathways_files,1)
        if any(contains(pathways_files(pathway_i).name,['SE_',slider_names{slider_i}]))
            pos_and_neg(slider_i) = 1;
        end
    end
end

