function ea_add_discfiber(~,~,leadgroup,resultfig)
% small function to open ea_discfiberexplorer. lead group can be a path to
% either lead group analysis file or saved fiber explorer.

if ismember('M', who('-file',leadgroup)) % Add new discfiber analysis
    exp = ea_discfiberexplorer(leadgroup,resultfig);
    tractId = exp.tractset.ID;

    ht=getappdata(resultfig,'ht');
    uipushtool(ht, 'CData', ea_get_icn('discfiber'),...
        'TooltipString', ['Explore fiber filtering analysis ',tractId],...
        'Tag', ['Explore fiber filtering analysis ',tractId],...
        'ClickedCallback', {@ea_add_discfiber,[fileparts(leadgroup),filesep,'disctracts',filesep,tractId,'.fibfilt'],resultfig});

    Tags = arrayfun(@(tool) tool.Tag, ht.Children, 'Uni', 0);
    isDiscFiberTool = contains(Tags, 'Explore fiber filtering analysis');
    isAddDiscFiberTool = contains(Tags, 'Add fiber filtering analysis');
    if any(isDiscFiberTool(2:end))
        insertInd = find(isDiscFiberTool(2:end),1);
    else
        insertInd = find(isAddDiscFiberTool,1) -1;
    end
    ht.Children=ht.Children([2:insertInd,1,insertInd+1:end]);
else % Load existing discfiber analysis
    ea_discfiberexplorer(leadgroup, resultfig);
end
