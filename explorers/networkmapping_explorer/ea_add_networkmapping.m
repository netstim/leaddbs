function ea_add_networkmapping(~,~,leadgroup,resultfig)
% small function to open ea_networkmappingexplorer. lead group can be a path to
% either lead group analysis file or saved fiber explorer.

if ismember('M', who('-file',leadgroup)) % Add new network mapping analysis
    exp = ea_networkmappingexplorer(leadgroup,resultfig);
    networkmappingId = exp.networkmapping.ID;

    ht=getappdata(resultfig,'ht');
    uipushtool(ht, 'CData', ea_get_icn('connectivities'),...
        'TooltipString', ['Explore DBS network mapping analysis ',networkmappingId],...
        'Tag', ['Explore DBS network mapping analysis ',networkmappingId],...
        'ClickedCallback', {@ea_add_networkmapping,[fileparts(leadgroup),filesep,'networkmapping',filesep,networkmappingId,'.netmap'],resultfig});

    Tags = arrayfun(@(tool) tool.Tag, ht.Children, 'Uni', 0);
    isNetworkMappingTool = contains(Tags, 'Explore DBS network mapping analysis');
    isAddNetworkMappingTool = contains(Tags, 'Add DBS network mapping analysis');
    if any(isNetworkMappingTool(2:end))
        insertInd = find(isNetworkMappingTool(2:end),1);
    else
        insertInd = find(isAddNetworkMappingTool,1) -1;
    end
    ht.Children=ht.Children([2:insertInd,1,insertInd+1:end]);
else % Load existing network mapping analysis
    ea_networkmappingexplorer(leadgroup, resultfig);
end
