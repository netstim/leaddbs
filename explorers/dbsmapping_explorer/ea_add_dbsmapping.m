function ea_add_dbsmapping(~,~,leadgroup,resultfig)
% small function to open ea_unifiedmappingexplorer. lead group can be a path to
% either lead group analysis file or saved fiber explorer.

if ismember('M', who('-file',leadgroup)) % Add new network mapping analysis
    exp = ea_dbsmappingexplorer(leadgroup,resultfig);
    dbsmappingId = exp.AnalysisIDEditField.Value;

    ht=getappdata(resultfig,'ht');
    uipushtool(ht, 'CData', ea_get_icn('connectivities'),...
        'TooltipString', ['Explore DBS mapping analysis ',dbsmappingId],...
        'Tag', ['Explore DBS mapping analysis ',dbsmappingId],...
        'ClickedCallback', {@ea_add_dbsmapping,[fileparts(leadgroup),filesep,'dbsmapping',filesep,dbsmappingId,'.dbsmap'],resultfig});

    Tags = arrayfun(@(tool) tool.Tag, ht.Children, 'Uni', 0);
    isUnifiedMappingTool = contains(Tags, 'Explore DBS mapping analysis');
    isAddUnifiedMappingTool = contains(Tags, 'Add DBS mapping analysis');
    if any(isUnifiedMappingTool(2:end))
        insertInd = find(isUnifiedMappingTool(2:end),1);
    else
        insertInd = find(isAddUnifiedMappingTool,1) -1;
    end
    ht.Children=ht.Children([2:insertInd,1,insertInd+1:end]);
else % Load existing unified mapping analysis
    ea_dbsmappingexplorer(leadgroup, resultfig);
end
