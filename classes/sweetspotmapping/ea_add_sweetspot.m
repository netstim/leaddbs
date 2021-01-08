function ea_add_sweetspot(~,~,leadgroup,resultfig)
% small function to open ea_discfiberexplorer. lead group can be a path to
% either lead group analysis file or saved fiber explorer.

if ismember('M', who('-file',leadgroup)) % Add new discfiber analysis
    exp = ea_sweetspotexplorer(leadgroup,resultfig);
    sweetId = exp.sweetspot.ID;

    ht=getappdata(resultfig,'ht');
    uipushtool(ht, 'CData', ea_get_icn('sweetspot'),...
        'TooltipString', ['Explore sweetspot analysis ',sweetId],...
        'Tag', ['Explore sweetspot analysis ',sweetId],...
        'ClickedCallback', {@ea_add_sweetspot,[fileparts(leadgroup),filesep,'sweetspots',filesep,sweetId,'.sweetspot'],resultfig});

    Tags = arrayfun(@(tool) tool.Tag, ht.Children, 'Uni', 0);
    isSweetSpotTool = contains(Tags, 'Explore sweetspot analysis');
    isAddSweetSpotTool = contains(Tags, 'Add sweetspot analysis');
    if any(isSweetSpotTool(2:end))
        insertInd = find(isSweetSpotTool(2:end),1);
    else
        insertInd = find(isAddSweetSpotTool,1) -1;
    end
    ht.Children=ht.Children([2:insertInd,1,insertInd+1:end]);
else % Load existing discfiber analysis
    ea_sweetspotexplorer(leadgroup, resultfig);
end
