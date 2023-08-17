function [connName, subset] = ea_connLabel2connName(connLabel)
% Return connectome name [and subset] based on connLabel

connName = [];
subset = [];

connBaseFolder = ea_getconnectomebase;
datasetInfo = ea_regexpdir(connBaseFolder, '^dataset_info\.json$');
[~, connNames] = fileparts(fileparts(datasetInfo));

for c=1:length(datasetInfo)
    try
        json = loadjson(datasetInfo{c});
    catch
        ea_cprintf('CmdWinErrors', 'Error loading connectome description file: %s ...\n', datasetInfo{c});
        continue;
    end

    if ~isfield(json, 'subsets')
        if strcmp(connLabel, ea_conn2connid(connNames{c}))
            connName = connNames{c};
            break;
        end
    else
        if startsWith(connLabel, ea_conn2connid(connNames{c}))
            connName = connNames{c};
            for s=1:length(json.subsets)
                if endsWith(connLabel, regexprep(json.subsets{s}.name, '[\W_]', ''))
                    subset = json.subsets{s}.name;
                    break;
                end
            end
            break;
        end
    end
end
